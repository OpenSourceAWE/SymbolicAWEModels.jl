# Issue draft — NetworkDynamics.jl

**Title:** Deeper-than-2 dependency chains / feedforward vertices — is this in scope?

Submit *after* the `cse` issue (`nd_cse_issue.md`), which is concrete, measured and comes with a patch. This one is a design question and lands better second.

Deliberately does not lead with "allow feedforward": asked that way it reads as a redesign request and invites a "that's what `ff_to_constraint` is for" reply. Leading with the concrete depth-4 chain and closing on the small ask gives an easy yes to (2) even if the answer to (1) is no.

---

A design question rather than a request, and related to #149 and #383.

We model kite structures with ModelingToolkit components and use ND as the runtime. Our dependency chain is depth 4: body states → body pose → attached ride-point position → segment force → back into the body. The coreloop's fixed vertex-g → edge-g → vertex-f order expresses depth 2, and the intermediate stages are genuine feedthrough vertices (outputs depending on inputs), which `construction.jl` rejects outside leaf+loopback. `ff_to_constraint` works but converts them to algebraic states, which we can't afford in the solve.

We work around it by collapsing those stages into their neighbours, and by zero-padding every component to one uniform I/O width. That is most of the complexity in our ND layer.

Two questions:

1. Is generalising the schedule (a build-time topological layering instead of the fixed 5-step loop) something you'd want in ND? A layered schedule is still statically known, so KA can launch per layer; the real cost is heterogeneous I/O width making buffers ragged rather than uniform-stride. Understood if that's a bigger change to the data model than you want.

2. If not, would you consider making `generate_io_function` public API? It already does the hard part — compiling a single MTK system with explicit I/O and unbound inputs — and it's usable standalone, with no `Network`, graph or coreloop. It's useful well beyond ND, and it would let projects like ours keep using your codegen with a stability guarantee instead of reaching into the extension.

Happy to write up our case in more detail if useful.
