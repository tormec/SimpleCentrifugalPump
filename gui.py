#!/usr/bin/env python3
"""Draw GUI for dimensioning simple centrifugal pump."""

import tkinter as tk
import tkinter.ttk as ttk
import simple_centrifugal_pump as scp
from simple_centrifugal_pump import Pre_Values as scp_pre


class Gui(tk.Frame):
    """Graphical User Interface for the script."""

    def __init__(self, master):
        """Initialize the interface."""
        tk.Frame.__init__(self)

        # notebook
        n = ttk.Notebook(master)
        n.grid(row=0, column=0)

        self.f_tab1 = tk.Frame(n)
        self.f_tab1.grid(row=0, column=0)

        self.f_tab2 = tk.Frame(n)
        self.f_tab2.grid(row=0, column=1)

        self.f_tab3 = tk.Frame(n)
        self.f_tab3.grid(row=0, column=3)

        n.add(self.f_tab1, text="Feasability Study")
        n.add(self.f_tab2, text="Pump Shaft")
        n.add(self.f_tab3, text="Impeller")

        # input variables
        f_input = tk.Frame(self.f_tab1)
        f_input.grid(row=0, column=0)

        l_flow = tk.Label(f_input, text="Flow [m^3/s]")
        l_flow.grid(row=0, column=0)
        self.flow = tk.DoubleVar()
        e_flow = tk.Entry(f_input, textvariable=self.flow, width=10)
        e_flow.grid(row=0, column=1, padx=5)

        l_head = tk.Label(f_input, text="Head [m]")
        l_head.grid(row=0, column=2)
        self.head = tk.DoubleVar()
        e_head = tk.Entry(f_input, textvariable=self.head, width=10)
        e_head.grid(row=0, column=3, padx=5)

        b_calc = tk.Button(f_input, text="Calculate",
                           command=self.print_pre_values)
        b_calc.grid(row=0, column=4)

        # table results feasability study
        f_fs_results = tk.Frame(self.f_tab1)
        f_fs_results.grid(row=1, column=0)

        unit_col = ["", "[rpm]", "", "", "", "[m/s]", "[mm]", "[mm]", "",
                    "[m]", ""]
        name_col = ["c. poles", "n", "k", "\u03C8", "\u03C6", "u1", "d1",
                    "b1", "b1/d1", "NPSHr", "\u03B7 tot"]
        title_col = list(zip(unit_col, name_col))

        for i, v in enumerate(title_col):
            l_title = tk.Label(f_fs_results, text=v[0] + "\n" + v[1])
            l_title.grid(row=0, column=i+1)

        self.choice = tk.IntVar()
        self.entry_result = []
        for r in range(len(scp.CPOLES)):
            for c in range(len(title_col)):
                pos = str(r) + "," + str(c)  # name widget
                var = tk.StringVar()
                self.entry_result.append([pos, var])
                e_result = tk.Entry(f_fs_results, width=8, textvariable=var)
                e_result.grid(row=r+1, column=c+1)
                r_result = tk.Radiobutton(f_fs_results, val=r,
                                          variable=self.choice)
                r_result.grid(row=r+1, column=0)

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
    """Execute GUI."""
    root = tk.Tk()

    app = Gui(master=root)
    app.master.title('Simple Centrifugal Pump')
    app.mainloop()

if __name__ == '__main__':
    main()
