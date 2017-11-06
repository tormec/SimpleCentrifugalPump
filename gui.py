#!/usr/bin/env python3
"""Draw GUI for edit/view data of geometrical dimensions for radial pump."""

import tkinter as tk
import tkinter.ttk as ttk
import simple_centrifugal_pump


class Gui(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self)  # initialize base class's init()
        self.grid()

        n = ttk.Notebook(master)
        n.grid(sticky=tk.N+tk.E+tk.S+tk.W, row=0, column=0)

        # frames for notebook
        self.f_tab1 = tk.Frame(n)
        self.f_tab1.grid(row=0, column=0)

        self.f_tab2 = tk.Frame(n)
        self.f_tab2.grid(row=0, column=1)

        self.f_tab3 = tk.Frame(n)
        self.f_tab3.grid(row=0, column=3)

        n.add(self.f_tab1, text="Feasability Study")
        n.add(self.f_tab2, text="Pump Shaft")
        n.add(self.f_tab3, text="Impeller")

        # entries for input variables
        l_flow = tk.Label(self.f_tab1, text="Flow [m^3/s]")
        l_flow.grid(row=0, column=0)
        self.flow = tk.DoubleVar()
        e_flow = tk.Entry(self.f_tab1, textvariable=self.flow, width=10)
        e_flow.grid(row=0, column=1, padx=5)

        l_head = tk.Label(self.f_tab1, text="Head [m]")
        l_head.grid(row=0, column=2)
        self.head = tk.DoubleVar()
        e_head = tk.Entry(self.f_tab1, textvariable=self.head, width=10)
        e_head.grid(row=0, column=3, padx=5)

        b_calc = tk.Button(self.f_tab1, text="Calculate", command=self.calc)
        b_calc.grid(row=0, column=4)

    def calc(self):
        print(self.flow.get())
        print(self.head.get())


def main():
    root = tk.Tk()

    app = Gui(master=root)
    app.master.title('Dimensioning Radial Pump')
    app.mainloop()

if __name__ == '__main__':
    main()
