/*
 * Copyright (C) 2014 The Android Open Source Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <fstream>
#include <stdint.h>

#include "builder.h"
#include "code_generator.h"
#include "compilers.h"
#include "driver/compiler_driver.h"
#include "driver/dex_compilation_unit.h"
#include "graph_visualizer.h"
#include "nodes.h"
#include "register_allocator.h"
#include "ssa_phi_elimination.h"
#include "ssa_liveness_analysis.h"
#include "utils/arena_allocator.h"

namespace art {

/**
 * Used by the code generator, to allocate the code in a vector.
 */
class CodeVectorAllocator FINAL : public CodeAllocator {
 public:
  CodeVectorAllocator() { }

  virtual uint8_t* Allocate(size_t size) {
    size_ = size;
    memory_.resize(size);
    return &memory_[0];
  }

  size_t GetSize() const { return size_; }
  const std::vector<uint8_t>& GetMemory() const { return memory_; }

 private:
  std::vector<uint8_t> memory_;
  size_t size_;

  DISALLOW_COPY_AND_ASSIGN(CodeVectorAllocator);
};

/**
 * If set to true, generates a file suitable for the c1visualizer tool and IRHydra.
 */
static bool kIsVisualizerEnabled = false;

/**
 * Filter to apply to the visualizer. Methods whose name contain that filter will
 * be in the file.
 */
static const char* kStringFilter = "";

class PassInfo;

class PassInfoPrinter : public ValueObject {
 public:
  PassInfoPrinter(HGraph* graph,
                  const char* method_name,
                  const CodeGenerator& codegen,
                  std::ostream* visualizer_output,
                  CompilerDriver* compiler_driver)
      : method_name_(method_name),
        timing_logger_enabled_(compiler_driver->GetDumpPasses()),
        timing_logger_(method_name, true, true),
        visualizer_enabled_(!compiler_driver->GetDumpCfgFileName().empty()),
        visualizer_(visualizer_output, graph, codegen) {
    if (strstr(method_name, kStringFilter) == nullptr) {
      timing_logger_enabled_ = visualizer_enabled_ = false;
    }
    if (visualizer_enabled_) {
      visualizer_.PrintHeader(method_name_);
    }
  }

  ~PassInfoPrinter() {
    if (timing_logger_enabled_) {
      LOG(INFO) << "TIMINGS " << method_name_;
      LOG(INFO) << Dumpable<TimingLogger>(timing_logger_);
    }
  }

 private:
  void StartPass(const char* pass_name) {
    // Dump graph first, then start timer.
    if (visualizer_enabled_) {
      visualizer_.DumpGraph(pass_name, /* is_after_pass */ false);
    }
    if (timing_logger_enabled_) {
      timing_logger_.StartTiming(pass_name);
    }
  }

  void EndPass(const char* pass_name) {
    // Pause timer first, then dump graph.
    if (timing_logger_enabled_) {
      timing_logger_.EndTiming();
    }
    if (visualizer_enabled_) {
      visualizer_.DumpGraph(pass_name, /* is_after_pass */ true);
    }
  }

  const char* method_name_;

  bool timing_logger_enabled_;
  TimingLogger timing_logger_;

  bool visualizer_enabled_;
  HGraphVisualizer visualizer_;

  friend PassInfo;

  DISALLOW_COPY_AND_ASSIGN(PassInfoPrinter);
};

class PassInfo : public ValueObject {
 public:
  PassInfo(const char *pass_name, PassInfoPrinter* pass_info_printer)
      : pass_name_(pass_name),
        pass_info_printer_(pass_info_printer) {
    pass_info_printer_->StartPass(pass_name_);
  }

  ~PassInfo() {
    pass_info_printer_->EndPass(pass_name_);
  }

 private:
  const char* const pass_name_;
  PassInfoPrinter* const pass_info_printer_;
};

class OptimizingCompiler FINAL : public Compiler {
 public:
  explicit OptimizingCompiler(CompilerDriver* driver);
  ~OptimizingCompiler();

  bool CanCompileMethod(uint32_t method_idx, const DexFile& dex_file, CompilationUnit* cu) const
      OVERRIDE;

  CompiledMethod* Compile(const DexFile::CodeItem* code_item,
                          uint32_t access_flags,
                          InvokeType invoke_type,
                          uint16_t class_def_idx,
                          uint32_t method_idx,
                          jobject class_loader,
                          const DexFile& dex_file) const OVERRIDE;

  CompiledMethod* TryCompile(const DexFile::CodeItem* code_item,
                             uint32_t access_flags,
                             InvokeType invoke_type,
                             uint16_t class_def_idx,
                             uint32_t method_idx,
                             jobject class_loader,
                             const DexFile& dex_file) const;

  CompiledMethod* JniCompile(uint32_t access_flags,
                             uint32_t method_idx,
                             const DexFile& dex_file) const OVERRIDE {
    return ArtQuickJniCompileMethod(GetCompilerDriver(), access_flags, method_idx, dex_file);
  }

  uintptr_t GetEntryPointOf(mirror::ArtMethod* method) const OVERRIDE
      SHARED_LOCKS_REQUIRED(Locks::mutator_lock_) {
    return reinterpret_cast<uintptr_t>(method->GetEntryPointFromQuickCompiledCodePtrSize(
        InstructionSetPointerSize(GetCompilerDriver()->GetInstructionSet())));
  }

  void InitCompilationUnit(CompilationUnit& cu) const OVERRIDE;

  void Init() OVERRIDE;

  void UnInit() const OVERRIDE;

  void MaybeRecordStat(MethodCompilationStat compilation_stat) const {
    if (compilation_stats_.get() != nullptr) {
      compilation_stats_->RecordStat(compilation_stat);
    }
  }

 private:
  // Whether we should run any optimization or register allocation. If false, will
  // just run the code generation after the graph was built.
  const bool run_optimizations_;

  // Optimize and compile `graph`.
  CompiledMethod* CompileOptimized(HGraph* graph,
                                   CodeGenerator* codegen,
                                   CompilerDriver* driver,
                                   const DexFile& dex_file,
                                   const DexCompilationUnit& dex_compilation_unit,
                                   PassInfoPrinter* pass_info) const;

  // Just compile without doing optimizations.
  CompiledMethod* CompileBaseline(CodeGenerator* codegen,
                                  CompilerDriver* driver,
                                  const DexCompilationUnit& dex_compilation_unit) const;

  std::unique_ptr<OptimizingCompilerStats> compilation_stats_;

  std::unique_ptr<std::ostream> visualizer_output_;

  // Delegate to Quick in case the optimizing compiler cannot compile a method.
  std::unique_ptr<Compiler> delegate_;

  DISALLOW_COPY_AND_ASSIGN(OptimizingCompiler);
};

static const int kMaximumCompilationTimeBeforeWarning = 100; /* ms */

OptimizingCompiler::OptimizingCompiler(CompilerDriver* driver)
    : Compiler(driver, kMaximumCompilationTimeBeforeWarning),
      run_optimizations_(
          (driver->GetCompilerOptions().GetCompilerFilter() != CompilerOptions::kTime)
          && !driver->GetCompilerOptions().GetDebuggable()),
      delegate_(Create(driver, Compiler::Kind::kQuick)) {}

void OptimizingCompiler::Init() {
  delegate_->Init();
  // Enable C1visualizer output. Must be done in Init() because the compiler
  // driver is not fully initialized when passed to the compiler's constructor.
  CompilerDriver* driver = GetCompilerDriver();
  const std::string cfg_file_name = driver->GetDumpCfgFileName();
  if (!cfg_file_name.empty()) {
    CHECK_EQ(driver->GetThreadCount(), 1U)
      << "Graph visualizer requires the compiler to run single-threaded. "
      << "Invoke the compiler with '-j1'.";
    visualizer_output_.reset(new std::ofstream(cfg_file_name));
  }
  if (driver->GetDumpStats()) {
    compilation_stats_.reset(new OptimizingCompilerStats());
  }
}

void OptimizingCompiler::UnInit() const {
  delegate_->UnInit();
}

OptimizingCompiler::~OptimizingCompiler() {
  if (compilation_stats_.get() != nullptr) {
    compilation_stats_->Log();
  }
}

void OptimizingCompiler::InitCompilationUnit(CompilationUnit& cu) const {
  delegate_->InitCompilationUnit(cu);
}

bool OptimizingCompiler::CanCompileMethod(uint32_t method_idx ATTRIBUTE_UNUSED,
                                          const DexFile& dex_file ATTRIBUTE_UNUSED,
                                          CompilationUnit* cu ATTRIBUTE_UNUSED) const {
  return true;
}

static bool IsInstructionSetSupported(InstructionSet instruction_set) {
  return instruction_set == kArm64
      || (instruction_set == kThumb2 && !kArm32QuickCodeUseSoftFloat)
      || instruction_set == kX86
      || instruction_set == kX86_64;
}

static bool CanOptimize(const DexFile::CodeItem& code_item) {
  // TODO: We currently cannot optimize methods with try/catch.
  return code_item.tries_size_ == 0;
}

static void RunOptimizations(HOptimization* optimizations[],
                             size_t length,
                             PassInfoPrinter* pass_info_printer) {
  for (size_t i = 0; i < length; ++i) {
    HOptimization* optimization = optimizations[i];
    {
      PassInfo pass_info(optimization->GetPassName(), pass_info_printer);
      optimization->Run();
    }
    optimization->Check();
  }
}

static void RunOptimizations(HGraph* graph,
                             CompilerDriver* driver,
                             OptimizingCompilerStats* stats,
                             const DexFile& dex_file,
                             const DexCompilationUnit& dex_compilation_unit,
                             PassInfoPrinter* pass_info_printer,
                             StackHandleScopeCollection* handles) {
  HDeadCodeElimination dce1(graph, stats);
  HDeadCodeElimination dce2(graph, stats, "dead_code_elimination_final");
  HConstantFolding fold1(graph);
  InstructionSimplifier simplify1(graph, stats);
  HBooleanSimplifier boolean_not(graph);

  HInliner inliner(graph, dex_compilation_unit, dex_compilation_unit, driver, stats);

  HConstantFolding fold2(graph);
  SideEffectsAnalysis side_effects(graph);
  GVNOptimization gvn(graph, side_effects);
  LICM licm(graph, side_effects);
  BoundsCheckElimination bce(graph);
  ReferenceTypePropagation type_propagation(graph, dex_file, dex_compilation_unit, handles);
  InstructionSimplifier simplify2(graph, stats, "instruction_simplifier_after_types");

  IntrinsicsRecognizer intrinsics(graph, dex_compilation_unit.GetDexFile(), driver);

  HOptimization* optimizations[] = {
    &intrinsics,
    &dce1,
    &fold1,
    &simplify1,
    // BooleanSimplifier depends on the InstructionSimplifier removing redundant
    // suspend checks to recognize empty blocks.
    &boolean_not,
    &inliner,
    &fold2,
    &side_effects,
    &gvn,
    &licm,
    &bce,
    &type_propagation,
    &simplify2,
    &dce2,
  };

  RunOptimizations(optimizations, arraysize(optimizations), pass_info_printer);
}

// The stack map we generate must be 4-byte aligned on ARM. Since existing
// maps are generated alongside these stack maps, we must also align them.
static ArrayRef<const uint8_t> AlignVectorSize(std::vector<uint8_t>& vector) {
  size_t size = vector.size();
  size_t aligned_size = RoundUp(size, 4);
  for (; size < aligned_size; ++size) {
    vector.push_back(0);
  }
  return ArrayRef<const uint8_t>(vector);
}

static void AllocateRegisters(HGraph* graph,
                              CodeGenerator* codegen,
                              PassInfoPrinter* pass_info_printer) {
  PrepareForRegisterAllocation(graph).Run();
  SsaLivenessAnalysis liveness(graph, codegen);
  {
    PassInfo pass_info(SsaLivenessAnalysis::kLivenessPassName, pass_info_printer);
    liveness.Analyze();
  }
  {
    PassInfo pass_info(RegisterAllocator::kRegisterAllocatorPassName, pass_info_printer);
    RegisterAllocator(graph->GetArena(), codegen, liveness).AllocateRegisters();
  }
}

CompiledMethod* OptimizingCompiler::CompileOptimized(HGraph* graph,
                                                     CodeGenerator* codegen,
                                                     CompilerDriver* compiler_driver,
                                                     const DexFile& dex_file,
                                                     const DexCompilationUnit& dex_compilation_unit,
                                                     PassInfoPrinter* pass_info_printer) const {
  StackHandleScopeCollection handles(Thread::Current());
  RunOptimizations(graph, compiler_driver, compilation_stats_.get(),
                   dex_file, dex_compilation_unit, pass_info_printer, &handles);

  AllocateRegisters(graph, codegen, pass_info_printer);

  CodeVectorAllocator allocator;
  codegen->CompileOptimized(&allocator);

  DefaultSrcMap src_mapping_table;
  if (compiler_driver->GetCompilerOptions().GetIncludeDebugSymbols()) {
    codegen->BuildSourceMap(&src_mapping_table);
  }
}

CompiledMethod* OptimizingCompiler::TryCompile(const DexFile::CodeItem* code_item,
                                               uint32_t access_flags,
                                               InvokeType invoke_type,
                                               uint16_t class_def_idx,
                                               uint32_t method_idx,
                                               jobject class_loader,
                                               const DexFile& dex_file) const {
  InstructionSet instruction_set = GetCompilerDriver()->GetInstructionSet();
  // Always use the thumb2 assembler: some runtime functionality (like implicit stack
  // overflow checks) assume thumb2.
  if (instruction_set == kArm) {
    instruction_set = kThumb2;
  }

  // Do not attempt to compile on architectures we do not support.
  if (instruction_set != kX86 && instruction_set != kX86_64 && instruction_set != kThumb2) {
    return nullptr;
  }

  DexCompilationUnit dex_compilation_unit(
    nullptr, class_loader, art::Runtime::Current()->GetClassLinker(), dex_file, code_item,
    class_def_idx, method_idx, access_flags,
    GetCompilerDriver()->GetVerifiedMethod(&dex_file, method_idx));

  // For testing purposes, we put a special marker on method names that should be compiled
  // with this compiler. This makes sure we're not regressing.
  bool shouldCompile = dex_compilation_unit.GetSymbol().find("00024opt_00024") != std::string::npos;
  bool shouldOptimize =
      dex_compilation_unit.GetSymbol().find("00024reg_00024") != std::string::npos;

  ArenaPool pool;
  ArenaAllocator arena(&pool);
  HGraphBuilder builder(&arena, &dex_compilation_unit, &dex_file, GetCompilerDriver());

  HGraph* graph = builder.BuildGraph(*code_item);
  if (graph == nullptr) {
    if (shouldCompile) {
      LOG(FATAL) << "Could not build graph in optimizing compiler";
    }
    return nullptr;
  }

  CodeGenerator* codegen = CodeGenerator::Create(&arena, graph, instruction_set);
  if (codegen == nullptr) {
    if (shouldCompile) {
      LOG(FATAL) << "Could not find code generator for optimizing compiler";
    }
    return nullptr;
  }

  HGraphVisualizer visualizer(
      visualizer_output_.get(), graph, kStringFilter, *codegen, dex_compilation_unit);
  visualizer.DumpGraph("builder");

  CodeVectorAllocator allocator;

  if (RegisterAllocator::CanAllocateRegistersFor(*graph, instruction_set)) {
    graph->BuildDominatorTree();
    graph->TransformToSSA();
    visualizer.DumpGraph("ssa");
    graph->FindNaturalLoops();

    SsaRedundantPhiElimination(graph).Run();
    SsaDeadPhiElimination(graph).Run();

    SsaLivenessAnalysis liveness(*graph, codegen);
    liveness.Analyze();
    visualizer.DumpGraph(kLivenessPassName);

    RegisterAllocator register_allocator(graph->GetArena(), codegen, liveness);
    register_allocator.AllocateRegisters();

    visualizer.DumpGraph(kRegisterAllocatorPassName);
    codegen->CompileOptimized(&allocator);
  } else if (shouldOptimize && RegisterAllocator::Supports(instruction_set)) {
    LOG(FATAL) << "Could not allocate registers in optimizing compiler";
  } else {
    codegen->CompileBaseline(&allocator);

    // Run these phases to get some test coverage.
    graph->BuildDominatorTree();
    graph->TransformToSSA();
    visualizer.DumpGraph("ssa");
    graph->FindNaturalLoops();
    SsaLivenessAnalysis liveness(*graph, codegen);
    liveness.Analyze();
    visualizer.DumpGraph(kLivenessPassName);
  }

  std::vector<uint8_t> mapping_table;
  codegen->BuildMappingTable(&mapping_table);
  std::vector<uint8_t> vmap_table;
  codegen->BuildVMapTable(&vmap_table);
  std::vector<uint8_t> gc_map;
  codegen->BuildNativeGCMap(&gc_map, dex_compilation_unit);

  return CompiledMethod::SwapAllocCompiledMethod(GetCompilerDriver(),
                                                 instruction_set,
                                                 ArrayRef<const uint8_t>(allocator.GetMemory()),
                                                 codegen->GetFrameSize(),
                                                 codegen->GetCoreSpillMask(),
                                                 0, /* FPR spill mask, unused */
                                                 ArrayRef<const uint8_t>(mapping_table),
                                                 ArrayRef<const uint8_t>(vmap_table),
                                                 ArrayRef<const uint8_t>(gc_map),
                                                 ArrayRef<const uint8_t>());
}

}  // namespace art
