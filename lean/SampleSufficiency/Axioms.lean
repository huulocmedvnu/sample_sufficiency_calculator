import SampleSufficiency.AngularError
import SampleSufficiency.Jacobian
-- Confirm the headline results depend only on standard Mathlib axioms
-- (propext, Classical.choice, Quot.sound) and nothing else.
#print axioms SampleSufficiency.trace_proj_conj
#print axioms SampleSufficiency.trace_proj_iso
#print axioms SampleSufficiency.trace_proj_diag
#print axioms SampleSufficiency.hasFDerivAt_normalize
