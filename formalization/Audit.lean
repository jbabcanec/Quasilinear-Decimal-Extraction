import BBP
import Lean.Util.CollectAxioms

/- Audit every public BBP declaration, including transitive dependencies.
The whitelist is the standard logical foundation used by mathlib; in
particular sorryAx and Lean.ofReduceBool are NOT permitted. -/

open Lean Elab Command in
run_cmd do
  let env ← getEnv
  let allowed : Array Name := #[``propext, ``Classical.choice, ``Quot.sound]
  let mut declarations : Array Name := #[]
  let mut theoremCount := 0
  for (name, info) in env.constants.toList do
    if (`BBP).isPrefixOf name then
      declarations := declarations.push name
      if info matches .thmInfo _ then
        theoremCount := theoremCount + 1
  declarations := declarations.qsort Name.quickLt
  for name in declarations do
    let axioms ← Lean.collectAxioms name
    for axiomName in axioms do
      unless allowed.contains axiomName do
        throwError "Unapproved axiom {axiomName} in {name}"
  unless theoremCount > 0 do
    throwError "No BBP theorems found; audit cannot pass an empty project"
  logInfo m!"AXIOM AUDIT PASS: {theoremCount} theorems; {declarations.size} declarations. Allowed: propext, Classical.choice, Quot.sound."
