"""
Transmission Line Voltage Regulation Calculator
-----------------------------------------------
A menu-driven program for engineering students to calculate the voltage
regulation and efficiency of a three-phase transmission line using the
short, medium (nominal-pi and nominal-T) and long line models.

All calculations are done per phase using ABCD (transmission) constants:
    Vs = A * Vr + B * Ir
    Is = C * Vr + D * Ir

Line constants (Z = total series impedance, Y = total shunt admittance):
    Short line        A = D = 1            B = Z                    C = 0
    Nominal-pi        A = D = 1 + ZY/2     B = Z                    C = Y(1 + ZY/4)
    Nominal-T         A = D = 1 + ZY/2     B = Z(1 + ZY/4)          C = Y
    Long line         A = D = cosh(gl)     B = Zc * sinh(gl)        C = sinh(gl) / Zc
        where  gamma = sqrt(z*y),  Zc = sqrt(z/y),  l = length

Voltage regulation:
    Vr(no load) = |Vs| / |A|
    % VR = (|Vr(no load)| - |Vr(full load)|) / |Vr(full load)| * 100

Transmission efficiency:
    Efficiency = Receiving-end power / Sending-end power * 100

Typical use:  Short line < 80 km,  Medium line 80 - 250 km,  Long line > 250 km.
A negative regulation (voltage rise) can occur with leading loads (Ferranti effect).
"""

import cmath
import math

SQRT3 = math.sqrt(3)

MODELS = {
    "short": "Short line",
    "pi": "Medium line (nominal-pi)",
    "t": "Medium line (nominal-T)",
    "long": "Long line (distributed)",
}


# ---------------------------------------------------------------- input helpers
def get_positive_float(prompt):
    """Keep asking until the user enters a valid positive number."""
    while True:
        try:
            value = float(input(prompt))
            if value <= 0:
                print("  Please enter a value greater than zero.")
                continue
            return value
        except ValueError:
            print("  Invalid input. Please enter a number.")


def get_nonneg_float(prompt):
    """Ask for a number that may be zero but not negative."""
    while True:
        try:
            value = float(input(prompt))
            if value < 0:
                print("  Please enter zero or a positive value.")
                continue
            return value
        except ValueError:
            print("  Invalid input. Please enter a number.")


def get_power_factor():
    """Ask for a power factor between 0 (exclusive) and 1 (inclusive)."""
    while True:
        pf = get_positive_float("Load power factor (0 to 1): ")
        if pf <= 1:
            return pf
        print("  Power factor cannot be greater than 1.")


def get_leading():
    """Ask whether the load power factor is lagging or leading."""
    while True:
        t = input("Power factor type - (L)agging or (E) leading: ").strip().lower()
        if t in ("l", "lagging"):
            return False
        if t in ("e", "leading"):
            return True
        print("  Please enter L or E.")


def get_line_data(need_shunt=False):
    """Ask for the line parameters (per phase, per km)."""
    length = get_positive_float("Line length (km): ")
    r = get_positive_float("Resistance r (ohm/km): ")
    x = get_positive_float("Inductive reactance x (ohm/km): ")
    if need_shunt:
        b = get_positive_float("Shunt susceptance b (microsiemens/km): ")
    else:
        b = get_nonneg_float("Shunt susceptance b (microsiemens/km, 0 if ignored): ")
    return length, r, x, b


def get_load_data():
    """Ask for the receiving-end load."""
    v_kv = get_positive_float("Receiving-end line voltage (kV): ")
    p_mw = get_positive_float("Receiving-end load power (MW, three-phase): ")
    pf = get_power_factor()
    leading = get_leading() if pf < 1 else False
    return v_kv, p_mw, pf, leading


# ------------------------------------------------------------------ core maths
def abcd_constants(model, length, r, x, b_us):
    """Return the complex A, B, C, D constants for the chosen line model."""
    z_per_km = complex(r, x)
    y_per_km = complex(0, b_us * 1e-6)
    z_total = z_per_km * length
    y_total = y_per_km * length

    if model == "short":
        return 1 + 0j, z_total, 0j, 1 + 0j

    if model == "pi":
        a = 1 + z_total * y_total / 2
        return a, z_total, y_total * (1 + z_total * y_total / 4), a

    if model == "t":
        a = 1 + z_total * y_total / 2
        return a, z_total * (1 + z_total * y_total / 4), y_total, a

    # long line (requires b > 0)
    gamma = cmath.sqrt(z_per_km * y_per_km)
    zc = cmath.sqrt(z_per_km / y_per_km)
    gl = gamma * length
    a = cmath.cosh(gl)
    return a, zc * cmath.sinh(gl), cmath.sinh(gl) / zc, a


def solve_line(abcd, v_kv, p_mw, pf, leading):
    """Solve the line with the receiving-end data. Returns a dictionary."""
    a, b, c, d = abcd

    vr = complex(v_kv * 1000 / SQRT3, 0)                    # phase voltage (reference)
    ir_mag = p_mw * 1e6 / (SQRT3 * v_kv * 1000 * pf)         # line current
    angle = math.acos(pf) * (1 if leading else -1)
    ir = cmath.rect(ir_mag, angle)

    vs = a * vr + b * ir
    i_s = c * vr + d * ir

    vr_no_load = abs(vs) / abs(a)
    regulation = (vr_no_load - abs(vr)) / abs(vr) * 100

    p_send = 3 * (vs * i_s.conjugate()).real
    p_recv = p_mw * 1e6
    efficiency = p_recv / p_send * 100

    return {
        "vs_line_kv": SQRT3 * abs(vs) / 1000,
        "vs_angle": math.degrees(cmath.phase(vs)),
        "is_mag": abs(i_s),
        "ir_mag": ir_mag,
        "pf_send": math.cos(cmath.phase(vs) - cmath.phase(i_s)),
        "vr_nl_line_kv": SQRT3 * vr_no_load / 1000,
        "regulation": regulation,
        "p_send_mw": p_send / 1e6,
        "loss_mw": (p_send - p_recv) / 1e6,
        "efficiency": efficiency,
    }


# ------------------------------------------------------------------- display
def show_abcd(abcd):
    """Print the ABCD constants in polar form."""
    for name, value in zip("ABCD", abcd):
        mag, ang = cmath.polar(value)
        if name == "C":
            print(f"  {name} = {mag:.6f} < {math.degrees(ang):.2f} deg S")
        elif name == "B":
            print(f"  {name} = {mag:.4f} < {math.degrees(ang):.2f} deg ohm")
        else:
            print(f"  {name} = {mag:.4f} < {math.degrees(ang):.2f} deg")


def show_results(title, abcd, result, v_kv):
    """Print full results for one model."""
    print(f"\n  ----- {title} -----")
    show_abcd(abcd)
    print(f"\n  Receiving-end line current = {result['ir_mag']:.2f} A")
    print(f"  Sending-end line voltage   = {result['vs_line_kv']:.3f} kV  (angle {result['vs_angle']:.2f} deg)")
    print(f"  Sending-end current        = {result['is_mag']:.2f} A")
    print(f"  Sending-end power factor   = {result['pf_send']:.4f}")
    print(f"  No-load receiving voltage  = {result['vr_nl_line_kv']:.3f} kV")
    print(f"  Full-load receiving voltage = {v_kv:.3f} kV")
    print(f"  VOLTAGE REGULATION         = {result['regulation']:.3f} %")
    print(f"  Sending-end power          = {result['p_send_mw']:.3f} MW")
    print(f"  Line power loss            = {result['loss_mw']:.3f} MW")
    print(f"  Transmission efficiency    = {result['efficiency']:.2f} %")
    if result["regulation"] < 0:
        print("  Note: negative regulation means the voltage RISES from no load to full load.")


def single_model(model):
    length, r, x, b = get_line_data(need_shunt=(model == "long"))
    v_kv, p_mw, pf, leading = get_load_data()
    abcd = abcd_constants(model, length, r, x, b)
    result = solve_line(abcd, v_kv, p_mw, pf, leading)
    show_results(MODELS[model], abcd, result, v_kv)


def compare_models():
    length, r, x, b = get_line_data(need_shunt=True)
    v_kv, p_mw, pf, leading = get_load_data()

    print("\n  ----- Model Comparison -----")
    print(f"  {'Model':28}{'VR (%)':>10}{'Vs (kV)':>11}{'Eff (%)':>10}")
    for key in ("short", "pi", "t", "long"):
        abcd = abcd_constants(key, length, r, x, b)
        res = solve_line(abcd, v_kv, p_mw, pf, leading)
        print(f"  {MODELS[key]:28}{res['regulation']:>10.3f}{res['vs_line_kv']:>11.3f}{res['efficiency']:>10.2f}")
    print("  The long-line model is the most accurate; short and medium models are approximations.")


def regulation_vs_power_factor():
    length, r, x, b = get_line_data(need_shunt=True)
    v_kv = get_positive_float("Receiving-end line voltage (kV): ")
    p_mw = get_positive_float("Receiving-end load power (MW, three-phase): ")

    model = "short" if length < 80 else "pi" if length <= 250 else "long"
    abcd = abcd_constants(model, length, r, x, b)

    print(f"\n  ----- Regulation vs Power Factor ({MODELS[model]}) -----")
    print(f"  {'Power factor':>14}{'Type':>10}{'VR (%)':>10}{'Vs (kV)':>11}")
    cases = [(0.7, False), (0.8, False), (0.9, False), (1.0, False), (0.9, True), (0.8, True), (0.7, True)]
    for pf, leading in cases:
        res = solve_line(abcd, v_kv, p_mw, pf, leading)
        kind = "unity" if pf == 1.0 else "leading" if leading else "lagging"
        print(f"  {pf:>14.2f}{kind:>10}{res['regulation']:>10.3f}{res['vs_line_kv']:>11.3f}")
    print("  Lagging loads give high regulation; leading loads can make it negative (Ferranti effect).")


def menu():
    print("\n" + "=" * 58)
    print("   TRANSMISSION LINE VOLTAGE REGULATION CALCULATOR")
    print("=" * 58)
    print(" 1. Short line (up to 80 km)")
    print(" 2. Medium line - nominal-pi (80 to 250 km)")
    print(" 3. Medium line - nominal-T (80 to 250 km)")
    print(" 4. Long line (above 250 km, distributed parameters)")
    print(" 5. Compare all models for the same line")
    print(" 6. Regulation vs load power factor table")
    print(" 0. Exit")
    print("-" * 58)


def main():
    while True:
        menu()
        choice = input("Enter your choice: ").strip()

        if choice == "1":
            single_model("short")
        elif choice == "2":
            single_model("pi")
        elif choice == "3":
            single_model("t")
        elif choice == "4":
            single_model("long")
        elif choice == "5":
            compare_models()
        elif choice == "6":
            regulation_vs_power_factor()
        elif choice == "0":
            print("\nThank you for using the calculator. Goodbye!")
            break
        else:
            print("  Invalid choice. Please select from the menu.")


if __name__ == "__main__":
    main()
