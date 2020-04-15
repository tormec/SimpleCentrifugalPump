#!/usr/bin/env python3
"""Draw GUI for dimensioning simple centrifugal pump."""

import tkinter as tk


class Gui(tk.Frame):
    """GUI for the script."""

    def __init__(self, master):
        """Initialize the interface."""
        tk.Frame.__init__(self)

        # input variables
        frm_input = tk.Frame(master)
        frm_input.grid(row=0, column=0)

        lbl_flow = tk.Label(frm_input, text="Flow [m\u00B3/s]", anchor="e")
        lbl_flow.grid(row=0, column=0, sticky="e, w")
        self.flow = tk.DoubleVar()
        ent_flow = tk.Entry(frm_input, textvariable=self.flow, width=10)
        ent_flow.grid(row=0, column=1, padx=5)

        lbl_head = tk.Label(frm_input, text="Head [m]", anchor="e")
        lbl_head.grid(row=0, column=2, sticky="e, w")
        self.head = tk.DoubleVar()
        ent_head = tk.Entry(frm_input, textvariable=self.head, width=10)
        ent_head.grid(row=0, column=3, padx=5)

        btn_calc = tk.Button(frm_input, text="Calculate",
                             command=self.print_results)
        btn_calc.grid(row=0, column=4)

        unit_col = ["", "[rpm]", "", "", "", "[m/s]", "[mm]", "[mm]", "",
                    "[m]", ""]
        name_col = ["c. poles", "n", "k", "\u03C8", "\u03C6", "u1", "d1",
                    "b1", "b1/d1", "NPSHr", "\u03B7 tot"]
        title_col = list(zip(unit_col, name_col))

        for i, v in enumerate(title_col):
            lbl_title = tk.Label(frm_input, text=v[0] + "\n" + v[1])
            lbl_title.grid(row=1, column=i)

        # self.choice = tk.IntVar()
        # self.entry_result = []
        # for r in range(len(scp.CPOLES)):
        #     for c in range(len(title_col)):
        #         pos = str(r) + "," + str(c)  # name widget
        #         var = tk.StringVar()
        #         self.entry_result.append([pos, var])
        #         ent_result = tk.Entry(frm_input, width=8,
        #                               textvariable=var)
        #         ent_result.grid(row=r+1, column=c+1)
        #         rdb_result = tk.Radiobutton(frm_input, val=r,
        #                                     variable=self.choice)
        #         rdb_result.grid(row=r+1, column=0)

    def print_results(self):
        pass


def main():
    """Run the GUI."""
    root = tk.Tk()

    app = Gui(master=root)
    app.master.title('Simple Centrifugal Pump')
    app.mainloop()


if __name__ == '__main__':
    main()
