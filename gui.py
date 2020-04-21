#!/usr/bin/env python3
"""Draw GUI for dimensioning simple centrifugal pump."""

import tkinter as tk
import simple_centrifugal_pump as scp


class Gui(tk.Frame):
    """GUI for the script."""

    def __init__(self, master):
        """Initialize the interface."""
        tk.Frame.__init__(self)

        # input
        frm_input = tk.Frame(master)
        frm_input.grid(row=0, column=0)

        # entry widget where input flow rate value
        lbl_flow = tk.Label(frm_input, text="flow rate [m^3/s]", anchor="e")
        lbl_flow.grid(row=0, column=0, sticky="ew")
        self.flow = tk.DoubleVar()
        ent_flow = tk.Entry(frm_input, textvariable=self.flow, width=10)
        ent_flow.grid(row=0, column=1, padx=5)

        # entry widget where input head value
        lbl_head = tk.Label(frm_input, text="head [m]", anchor="e")
        lbl_head.grid(row=0, column=2, sticky="ew")
        self.head = tk.DoubleVar()
        ent_head = tk.Entry(frm_input, textvariable=self.head, width=10)
        ent_head.grid(row=0, column=3, padx=5)

        # button widget which execute the main script
        btn_calc = tk.Button(frm_input, text="calculate",
                             command=self.print_results)
        btn_calc.grid(row=0, column=4)

        # output
        frm_output = tk.Frame(master)
        frm_output.grid(row=1, column=0)

        # text widget where to print results
        self.txt_results = tk.Text(frm_output, state="disabled")
        self.txt_results.grid(row=1, column=0)
        # add vertical scrollbar
        vbar = tk.Scrollbar(frm_output, command=self.txt_results.yview)
        vbar.grid(row=1, column=1, sticky="ns")
        self.txt_results["yscroll"] = vbar.set
        # allow to copy the selection with ctrl+c
        self.txt_results.bind("<1>",
                              lambda event: self.txt_results.focus_set())

    def print_results(self):
        """Print the results."""
        results = scp.main(flow=self.flow.get(), head=self.head.get())

        self.txt_results["state"] = "normal"
        self.txt_results.insert("end", results)
        self.txt_results["state"] = "disabled"


def main():
    """Run the GUI."""
    root = tk.Tk()

    app = Gui(master=root)
    app.master.title('Simple Centrifugal Pump')
    app.mainloop()


if __name__ == '__main__':
    main()
