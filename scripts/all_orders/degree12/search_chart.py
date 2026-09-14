#!/usr/bin/env python3
import json
from pathlib import Path
from chiral4form.degree12_chart_search import search_charts

if __name__=="__main__":
    r=search_charts()
    print("\n===== DEGREE-12 CHART SEARCH =====")
    print("status:",r["status"])
    print("attempts:",r["attempt_count"])
    if r["status"]=="exact_chart_found":
        print("selected free columns:",r["selected_free"])
        print("PASS: exact degree-12 QQ kernel chart reconstructed")
    else:
        print("best free columns:",r.get("best_free"))
        print("best unresolved:",r.get("best_unresolved"))
        print("NEXT:",r["next_action"])
