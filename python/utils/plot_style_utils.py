import matplotlib as plt
import ROOT

def draw_title_text(x_min, x_max, y_min, y_max, text):
    pt = ROOT.TPaveText(x_min, y_min, x_max, y_max, "brNDC")
    pt.SetTextFont(42)
    pt.AddText(f"#scale[1.35]{{{text}}}")
    pt.SetFillStyle(0)   # Transparent background
    pt.SetBorderSize(0)  # No border
    pt.Draw("same")

def draw_collision_energy(x_min, x_max, y_min, y_max, system, energy):
    pt = ROOT.TPaveText(x_min, y_min, x_max, y_max, "brNDC")
    pt.SetTextFont(42)

    if system == 'pp':
        pt.AddText(f"#scale[1.2]{{{system}}}, #sqrt{{#it{{s}}}} = {energy} TeV")
    else:
        pt.AddText(f"#scale[1.2]{{{system}}}, #sqrt{{#it{{s}}_{{NN}}}} = {energy} TeV")

    pt.SetFillStyle(0)   # Transparent background
    pt.SetBorderSize(0)  # No border
    pt.Draw("same")

def draw_d_meson_decay(x_min, x_max, y_min, y_max, meson):
    pt = ROOT.TPaveText(x_min, y_min, x_max, y_max, "brNDC")
    pt.SetTextFont(42)

    if meson == 'D0':
        pt.AddText("D^{0} #rightarrow K^{#minus}#pi^{+} and charge conj.")
    elif meson == 'Dplus':
        pt.AddText("D^{+} #rightarrow K^{#minus}#pi^{+}#pi^{+} and charge conj.")
    #TODO: add more mesons

    pt.SetFillStyle(0)   # Transparent background
    pt.SetBorderSize(0)  # No border
    pt.Draw("same")

def draw_title_collision_energy(x_min, x_max, y_min, y_max, title, system, energy):
    draw_title_text(x_min, x_max, y_min, y_max, title)
    draw_collision_energy(x_min, x_max, y_min-0.1, y_max-0.1, system, energy)


def draw_all_info(x_min, x_max, y_min, y_max, title, system, energy, meson):
    print('Drawing all info')
    draw_title_collision_energy(x_min, x_max, y_min, y_max, title, system, energy)
    draw_d_meson_decay(x_min, x_max, y_min-0.15, y_max-0.15, meson)