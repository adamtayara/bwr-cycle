"""
BWR power-cycle model (direct-cycle boiling water reactor).

Compares two cycles using IAPWS-IF97 steam tables:
  1. Simple Rankine: reactor steam expands straight to the condenser.
  2. BWR-style cycle: HP turbine -> moisture separator -> LP turbine,
     with the separator drains and an HP-exhaust extraction feeding an
     open feedwater heater.

Sweeps reactor (steam dome) pressure and crossover pressure, reports
thermal efficiency and turbine-exhaust moisture, and scales the result
to an 870 MWt reactor (roughly BWRX-300 thermal rating).

Run:  python bwr_cycle.py
"""

import numpy as np
import matplotlib.pyplot as plt
from iapws import IAPWS97

# ---------------- inputs ----------------
P_REACTOR = 7.07       # MPa, typical BWR steam dome pressure
P_COND = 0.005         # MPa (5 kPa condenser)
ETA_TURB = 0.85        # isentropic efficiency, HP and LP turbines
ETA_PUMP = 0.80        # isentropic efficiency, pumps
ETA_GEN = 0.98         # generator efficiency
Q_REACTOR_MW = 870.0   # reactor thermal power, MWt
MIN_EXHAUST_X = 0.88   # rough lower limit on LP exhaust quality (blade erosion)


def sat_vap(P):
    return IAPWS97(P=P, x=1.0)


def sat_liq(P):
    return IAPWS97(P=P, x=0.0)


def expand(state_in, P_out, eta):
    """Turbine expansion with isentropic efficiency. Returns (h_out, x_out)."""
    h_s = IAPWS97(P=P_out, s=state_in.s).h
    h_out = state_in.h - eta * (state_in.h - h_s)
    out = IAPWS97(P=P_out, h=h_out)
    return out.h, out.x


def pump(liq_in, P_out, eta):
    """Pump work with isentropic efficiency. Returns h_out."""
    h_s = IAPWS97(P=P_out, s=liq_in.s).h
    return liq_in.h + (h_s - liq_in.h) / eta


def simple_rankine(Pr, Pc=P_COND):
    s1 = sat_vap(Pr)
    h2, x2 = expand(s1, Pc, ETA_TURB)
    s3 = sat_liq(Pc)
    h4 = pump(s3, Pr, ETA_PUMP)
    w_net = (s1.h - h2) - (h4 - s3.h)
    q_in = s1.h - h4
    return {"eta": w_net / q_in, "x_exhaust": x2, "w_net": w_net, "q_in": q_in}


def bwr_cycle(Pr, Pm, Pc=P_COND):
    """Per kg of reactor steam. Pm = HP exhaust / crossover pressure."""
    s1 = sat_vap(Pr)                         # reactor outlet
    h2, x2 = expand(s1, Pm, ETA_TURB)        # HP turbine exhaust (wet)

    s3 = sat_vap(Pm)                         # separator vapor outlet
    hf_m = sat_liq(Pm).h                     # separator drain
    h4, x4 = expand(s3, Pc, ETA_TURB)        # LP turbine exhaust

    s5 = sat_liq(Pc)                         # condensate
    h6 = pump(s5, Pm, ETA_PUMP)              # condensate pump
    s7 = sat_liq(Pm)                         # open FWH outlet
    h8 = pump(s7, Pr, ETA_PUMP)              # feedwater pump

    # Open FWH energy balance (1 kg out):
    # f*h2 + (1-f)(1-x2)*hf_m + (1-f)*x2*h6 = h7  -> solve for extraction f
    a = (1 - x2) * hf_m + x2 * h6
    f = (s7.h - a) / (h2 - a)

    m_lp = (1 - f) * x2                      # kg to LP turbine per kg reactor steam
    w_hp = s1.h - h2
    w_lp = m_lp * (s3.h - h4)
    w_pumps = m_lp * (h6 - s5.h) + (h8 - s7.h)
    w_net = w_hp + w_lp - w_pumps
    q_in = s1.h - h8
    return {"eta": w_net / q_in, "x_exhaust": x4, "x_hp": x2,
            "extraction": f, "w_net": w_net, "q_in": q_in, "T_fw": IAPWS97(P=Pr, h=h8).T - 273.15}


def main():
    # 1) Baseline at nominal pressure
    simple = simple_rankine(P_REACTOR)

    # 2) Best crossover pressure at nominal reactor pressure
    Pm_grid = np.linspace(0.3, 3.0, 55)
    eta_vs_pm = [bwr_cycle(P_REACTOR, pm)["eta"] for pm in Pm_grid]
    Pm_best = Pm_grid[int(np.argmax(eta_vs_pm))]
    best = bwr_cycle(P_REACTOR, Pm_best)

    # 3) Reactor pressure sweep (at best crossover pressure)
    Pr_grid = np.linspace(5.0, 8.0, 31)
    eta_simple = [simple_rankine(p)["eta"] for p in Pr_grid]
    eta_bwr = [bwr_cycle(p, Pm_best)["eta"] for p in Pr_grid]
    x_simple = [simple_rankine(p)["x_exhaust"] for p in Pr_grid]
    x_bwr = [bwr_cycle(p, Pm_best)["x_exhaust"] for p in Pr_grid]

    # 4) Scale to reactor
    P_elec = Q_REACTOR_MW * best["eta"] * ETA_GEN
    m_dot = Q_REACTOR_MW * 1e3 / best["q_in"]   # kg/s (q_in in kJ/kg)

    print("=== Nominal case: P_reactor = %.2f MPa, P_cond = %.0f kPa ===" % (P_REACTOR, P_COND * 1e3))
    print("Simple Rankine : eta = %.1f%%, LP exhaust quality = %.3f" % (simple["eta"] * 100, simple["x_exhaust"]))
    print("BWR cycle      : eta = %.1f%%, LP exhaust quality = %.3f" % (best["eta"] * 100, best["x_exhaust"]))
    print("  best crossover pressure = %.2f MPa" % Pm_best)
    print("  HP exhaust quality      = %.3f" % best["x_hp"])
    print("  HP extraction fraction  = %.3f" % best["extraction"])
    print("  feedwater temperature   = %.0f C" % best["T_fw"])
    print("Efficiency gain vs simple: %+.1f points" % ((best["eta"] - simple["eta"]) * 100))
    print("Exhaust moisture: %.1f%% -> %.1f%%" % ((1 - simple["x_exhaust"]) * 100, (1 - best["x_exhaust"]) * 100))
    print("\n=== Scaled to %.0f MWt ===" % Q_REACTOR_MW)
    print("Electric output ~ %.0f MWe, steam flow ~ %.0f kg/s" % (P_elec, m_dot))
    print("\nReactor pressure sweep 5-8 MPa: BWR eta %.1f%% -> %.1f%%" % (eta_bwr[0] * 100, eta_bwr[-1] * 100))

    # Plots
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    ax[0].plot(Pm_grid, np.array(eta_vs_pm) * 100)
    ax[0].axvline(Pm_best, ls="--", c="gray")
    ax[0].set(xlabel="Crossover pressure (MPa)", ylabel="Thermal efficiency (%)",
              title="Efficiency vs crossover pressure")
    ax[1].plot(Pr_grid, np.array(eta_simple) * 100, label="Simple Rankine")
    ax[1].plot(Pr_grid, np.array(eta_bwr) * 100, label="BWR cycle (MS + FWH)")
    ax[1].set(xlabel="Reactor pressure (MPa)", ylabel="Thermal efficiency (%)",
              title="Efficiency vs reactor pressure")
    ax[1].legend()
    ax[2].plot(Pr_grid, x_simple, label="Simple Rankine")
    ax[2].plot(Pr_grid, x_bwr, label="BWR cycle (MS + FWH)")
    ax[2].axhline(MIN_EXHAUST_X, ls="--", c="red", label="~erosion limit")
    ax[2].set(xlabel="Reactor pressure (MPa)", ylabel="LP exhaust quality",
              title="Turbine exhaust moisture")
    ax[2].legend()
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig("bwr_cycle_results.png", dpi=150)
    print("\nSaved bwr_cycle_results.png")


if __name__ == "__main__":
    main()
