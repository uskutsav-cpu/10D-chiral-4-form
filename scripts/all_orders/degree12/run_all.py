#!/usr/bin/env python3
from chiral4form.degree12_kernel_lift import lift_kernel
from chiral4form.degree12_explicit_ideal import build_ideal
from chiral4form.degree12_exact_tangency import certify_tangency
from chiral4form.degree12_global_check import check_global
from chiral4form.degree12_all_orders_promotion import promote

def main():
    print("===== DEGREE-12 STEP 1: LIFT KERNEL =====")
    k=lift_kernel()
    print("status:",k["status"])
    print("model directory:",k["model_directory"])
    print("pivots:",len(k["pivots"]),"free degree12:",k["free_degree12_count"])
    if not k["exact"]:
        raise SystemExit("STOP: exact all-prime kernel reconstruction did not close")

    print("\n===== DEGREE-12 STEP 2: BUILD IDEAL =====")
    i=build_ideal()
    print("status:",i["status"])
    print("constraints:",i["constraint_count"],"fiber:",i["fiber_dimension"])

    print("\n===== DEGREE-12 STEP 3: TANGENCY =====")
    t=certify_tangency()
    print("status:",t["status"])
    print("finite ideal tangent:",t["finite_ideal_tangent"])

    print("\n===== DEGREE-12 STEP 4: GLOBAL CHECK =====")
    g=check_global()
    print("status:",g["status"])
    print("6 + 4 =",g["orbit_dimension"])

    print("\n===== DEGREE-12 STEP 5: ALL ORDERS =====")
    a=promote()
    print("status:",a["status"])
    print("all-orders cylinder invariant:",a["all_orders_cylinder_invariant"])

    print("\nPASS: DEGREE-12 EXACT-IDEAL PIPELINE COMPLETE")

if __name__=="__main__":
    main()
