#!/usr/bin/env python3
"""Draw GUI for dimensioning simple centrifugal pump."""

import tkinter as tk
import tkinter.ttk as ttk
import simple_centrifugal_pump as scp
from simple_centrifugal_pump import Pre_Values as scp_pre


class Gui(tk.Frame):
    """GUI for the script."""

    def __init__(self, master):
        """Initialize the interface."""
        tk.Frame.__init__(self)

        # notebook
        ntb = ttk.Notebook(master)
        ntb.grid(row=0, column=0)

        self.frm_tab1 = tk.Frame(ntb)
        self.frm_tab1.grid(row=0, column=0)

        self.frm_tab2 = tk.Frame(ntb)
        self.frm_tab2.grid(row=0, column=1)

        self.frm_tab3 = tk.Frame(ntb)
        self.frm_tab3.grid(row=0, column=3)

        ntb.add(self.frm_tab1, text="Feasability Study")
        ntb.add(self.frm_tab2, text="Pump Shaft")
        ntb.add(self.frm_tab3, text="Impeller")

        # input variables
        frm_input = tk.Frame(master)
        frm_input.grid(row=0, column=1)

        lbl_flow = tk.Label(frm_input, text="Flow [m^3/s]", anchor="e")
        lbl_flow.grid(row=0, column=0, sticky="e, w")
        self.flow = tk.DoubleVar()
        ent_flow = tk.Entry(frm_input, textvariable=self.flow, width=10)
        ent_flow.grid(row=0, column=1, padx=5)

        lbl_head = tk.Label(frm_input, text="Head [m]", anchor="e")
        lbl_head.grid(row=1, column=0, sticky="e, w")
        self.head = tk.DoubleVar()
        ent_head = tk.Entry(frm_input, textvariable=self.head, width=10)
        ent_head.grid(row=1, column=1, padx=5)

        lbl_hz = tk.Label(frm_input, text="f [Hz]", anchor="e")
        lbl_hz.grid(row=2, column=0, sticky="e, w")
        self.hz = tk.DoubleVar()
        self.hz.set("50")
        ent_hz = tk.Entry(frm_input, textvariable=self.hz, width=10)
        ent_hz.grid(row=2, column=1, padx=5)

        lbl_slip = tk.Label(frm_input, text="slip [%]", anchor="e")
        lbl_slip.grid(row=3, column=0, sticky="e, w")
        self.slip = tk.DoubleVar()
        self.slip.set("3")
        ent_slip = tk.Entry(frm_input, textvariable=self.slip, width=10)
        ent_slip.grid(row=3, column=1, padx=5)

        lbl_tau_adm = tk.Label(frm_input, text="\u03C4 adm [MPa]", anchor="e")
        lbl_tau_adm.grid(row=4, column=0, sticky="e, w")
        self.tau_adm = tk.DoubleVar()
        ent_tau_adm = tk.Entry(frm_input, textvariable=self.tau_adm, width=10)
        ent_tau_adm.grid(row=4, column=1, padx=5)

        lbl_thk = tk.Label(frm_input, text="blade thk. [m]", anchor="e")
        lbl_thk.grid(row=5, column=0, sticky="e, w")
        self.thk = tk.DoubleVar()
        self.thk.set(".003")
        ent_thk = tk.Entry(frm_input, textvariable=self.thk, width=10)
        ent_thk.grid(row=5, column=1, padx=5)

        lbl_lm = tk.Label(frm_input, text="lm", anchor="e")
        lbl_lm.grid(row=6, column=0, sticky="e, w")
        self.lm = tk.DoubleVar()
        self.lm.set(".04")
        ent_lm = tk.Entry(frm_input, textvariable=self.lm, width=10)
        ent_lm.grid(row=6, column=1, padx=5)

        lbl_lw = tk.Label(frm_input, text="lw", anchor="e")
        lbl_lw.grid(row=7, column=0, sticky="e, w")
        self.lw = tk.DoubleVar()
        self.lw.set(".50")
        ent_lw = tk.Entry(frm_input, textvariable=self.lw, width=10)
        ent_lw.grid(row=7, column=1, padx=5)

        lbl_km = tk.Label(frm_input, text="km", anchor="e")
        lbl_km.grid(row=8, column=0, sticky="e, w")
        self.km = tk.DoubleVar()
        self.km.set("1.2")
        ent_km = tk.Entry(frm_input, textvariable=self.km, width=10)
        ent_km.grid(row=8, column=1, padx=5)

        lbl_eta_vol = tk.Label(frm_input, text="\u03B7 vol", anchor="e")
        lbl_eta_vol.grid(row=9, column=0, sticky="e, w")
        self.eta_vol = tk.DoubleVar()
        ent_eta_vol = tk.Entry(frm_input, textvariable=self.eta_vol, width=10)
        ent_eta_vol.grid(row=9, column=1, padx=5)

        lbl_eta_idr = tk.Label(frm_input, text="\u03B7 idr", anchor="e")
        lbl_eta_idr.grid(row=10, column=0, sticky="e, w")
        self.eta_idr = tk.DoubleVar()
        ent_eta_idr = tk.Entry(frm_input, textvariable=self.eta_idr, width=10)
        ent_eta_idr.grid(row=10, column=1, padx=5)

        lbl_d2 = tk.Label(frm_input, text="d2 [m]", anchor="e")
        lbl_d2.grid(row=11, column=0, sticky="e, w")
        self.d2 = tk.DoubleVar()
        ent_d2 = tk.Entry(frm_input, textvariable=self.d2, width=10)
        ent_d2.grid(row=11, column=1, padx=5)

        lbl_gamma_2 = tk.Label(frm_input, text="\u03b3 [deg]", anchor="e")
        lbl_gamma_2.grid(row=12, column=0, sticky="e, w")
        self.gamma_2 = tk.IntVar()
        ent_gamma_2 = tk.Entry(frm_input, textvariable=self.gamma_2, width=10)
        ent_gamma_2.grid(row=12, column=1, padx=5)

        lbl_z = tk.Label(frm_input, text="number of blade", anchor="e")
        lbl_z.grid(row=13, column=0, sticky="e, w")
        self.z = tk.IntVar()
        ent_z = tk.Entry(frm_input, textvariable=self.z, width=10)
        ent_z.grid(row=13, column=1, padx=5)

        btn_calc = tk.Button(frm_input, text="Calculate",
                             command=self.print_pre_values)
        btn_calc.grid(row=14, column=0, columnspan=2)

        # table results feasability study
        frm_fs_results = tk.Frame(self.frm_tab1)
        frm_fs_results.grid(row=1, column=0)

        unit_col = ["", "[rpm]", "", "", "", "[m/s]", "[mm]", "[mm]", "",
                    "[m]", ""]
        name_col = ["c. poles", "n", "k", "\u03C8", "\u03C6", "u1", "d1",
                    "b1", "b1/d1", "NPSHr", "\u03B7 tot"]
        title_col = list(zip(unit_col, name_col))

        for i, v in enumerate(title_col):
            lbl_title = tk.Label(frm_fs_results, text=v[0] + "\n" + v[1])
            lbl_title.grid(row=0, column=i+1)

        self.choice = tk.IntVar()
        self.entry_result = []
        for r in range(len(scp.CPOLES)):
            for c in range(len(title_col)):
                pos = str(r) + "," + str(c)  # name widget
                var = tk.StringVar()
                self.entry_result.append([pos, var])
                ent_result = tk.Entry(frm_fs_results, width=8,
                                      textvariable=var)
                ent_result.grid(row=r+1, column=c+1)
                rdb_result = tk.Radiobutton(frm_fs_results, val=r,
                                            variable=self.choice)
                rdb_result.grid(row=r+1, column=0)

    def print_pre_values(self):

        fs_rpm = scp_pre.rotational_speed(self.slip, self.hz)
        fs_k_num = scp_pre.type_number(self.fs_rpm, self.flow, self.head)
        fs_u1 = scp_pre._circumferential_velocity_1(self.head, self.fs_psi)
        fs_d1 = scp_pre._diameter_1(self.fs_u1, self.fs_rpm)
        fs_b1 = scp_pre._width_1(self.fs_u1, self.fs_d1, self.flow,
                                 self.fs_phi)
        fs_bd1 = scp_pre.width_over_diameter_1(self.fs_b1, self.fs_d1)
        fs_npsh_r = scp_pre.npsh_r(self.fs_k_num, self.head)

        values = zip(self.fs_cpoles, self.fs_rpm, self.fs_k_num, self.fs_psi,
                     self.fs_phi, self.fs_u1, self.fs_d1, self.fs_b1,
                     self.fs_bd1, self.fs_npsh_r, self.fs_eta)
        values = list(values)
        for i in iter(self.entry_result):  # i = ["r,c", var]
            rc = i[0].split(",")  # rc = ["r", "c"]
            r = int(rc[0])  # row number
            c = int(rc[1])  # column number
            i[1].set(round(values[r][c], 5))


def main():
    """Run the GUI."""
    root = tk.Tk()

    app = Gui(master=root)
    app.master.title('Simple Centrifugal Pump')
    app.mainloop()

if __name__ == '__main__':
    main()
