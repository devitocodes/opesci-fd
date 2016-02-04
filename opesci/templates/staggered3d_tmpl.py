from regular3d_tmpl import *

class Staggered3DTemplate(Regular3DTemplate):
    # Order of names in the following list is important. The resulting code blocks would be placed in the same order as they appear here
    _template_methods = ['includes', 'grid_structure', 'convergence_structure', 'profiling_function', 'execute', 'convergence_function', 'freemem', 'main']
    __convergence_structure_name = 'OpesciConvergence'
    
    def convergence_structure(self):
        return cgen.Extern("C", cgen.Struct(self.__convergence_structure_name, [self.grid.define_convergence]))
    
    def execute_parallel_block(self):
        statements = []
        if self.profiling:
            if self.numevents_papi>0:
                statements += [self.grid.define_papi_events]
                statements.append(cgen.Statement("opesci_papi_start_counters(numevents, events)"))
            else:
                statements.append(cgen.Value("float", "real_time"))
                statements.append(cgen.Value("float", "proc_time"))
                statements.append(cgen.Value("float", "mflops"))
                statements.append(cgen.Value("long long", "flpins"))
                statements.append(cgen.Statement("opesci_flops(&real_time, &proc_time, &flpins, &mflops)"))
        statements.append(self.grid.initialise)
        statements.append(self.grid.initialise_bc)
        statements.append(self.execute_time_loop())

        if self.profiling:
            if self.numevents_papi>0:
                statements.append(cgen.Statement("opesci_papi_read_counters(numevents, counters)"))
                statements.append(cgen.Pragma("omp critical"))
                statements.append(cgen.Block(self.grid.sum_papi_events()))
            else:
                statements.append(cgen.Statement("opesci_flops(&real_time, &proc_time, &flpins, &mflops)"))
                statements.append(cgen.Pragma("omp critical"))
                critical_block = []
                critical_block.append(cgen.Assign("profiling->g_rtime", "fmax(profiling->g_rtime, real_time)"))
                critical_block.append(cgen.Assign("profiling->g_ptime", "fmax(profiling->g_ptime, proc_time)"))
                critical_block.append(cgen.Statement("profiling->g_mflops += mflops;"))
                statements.append(cgen.Block(critical_block))
        print statements
        return [cgen.Pragma("omp parallel"), cgen.Block(statements)]
    
    def execute_time_loop(self):
        statements = []
        statements.append(self.grid.time_stepping)
        if self.pluto:
            statements.append(cgen.Pragma("scop"))
        statements.append(self.grid.stress_loop)
        if self.pluto:
            statements.append(cgen.Pragma("endscop"))
        statements.append(self.grid.stress_bc)
        
        if self.pluto:
            statements.append(cgen.Pragma("scop"))
        statements.append(self.grid.velocity_loop)
        if self.pluto:
            statements.append(cgen.Pragma("endscop"))
        statements.append(self.grid.velocity_bc)
        output_step = self.grid.output_step
        if output_step:
            statements.append(output_step)
        result = cgen.For(cgen.InlineInitializer(cgen.Value("int", "_ti"), 0), "_ti < ntsteps", "_ti++", cgen.Block(statements))
        return result

    def convergence_function(self):
        statements = []
        statements.append(self.grid.define_constants)
        statements.append(self.grid.load_fields)
        statements.append(self.grid.converge_test)
        statements.append(cgen.Statement("return 0"))
        return cgen.FunctionBody(cgen.Extern("C", cgen.FunctionDeclaration(cgen.Value('int', 'opesci_convergence'), [cgen.Pointer(cgen.Value(self._grid_structure_name, "grid")), cgen.Pointer(cgen.Value(self.__convergence_structure_name, "conv"))])), cgen.Block(statements))
    
    def main(self):
        statements = []
        statements.append(cgen.Value(self._grid_structure_name, "grid"))
        statements.append(cgen.Value(self.__convergence_structure_name, "conv"))
        statements.append(cgen.Value(self._profiling_structure_name, "profiling"))
        statements.append(cgen.Statement("opesci_execute(&grid, &profiling)"))
        statements.append(cgen.Statement("opesci_convergence(&grid, &conv)"))
        statements.append(cgen.Statement("opesci_free(&grid)"))
        statements.append(self.grid.print_convergence)
        if self.profiling:
            statements.append(cgen.Statement('printf("PAPI:: Max real_time: %f (sec)\\n", profiling.g_rtime)'))
            statements.append(cgen.Statement('printf("PAPI:: Max proc_time: %f (sec)\\n", profiling.g_ptime)'))
            statements.append(cgen.Statement('printf("PAPI:: Total MFlops/s: %f\\n", profiling.g_mflops)'))
        statements.append(cgen.Statement('return 0'))
        return cgen.FunctionBody(cgen.FunctionDeclaration(cgen.Value('int', 'main'), []), cgen.Block(statements))