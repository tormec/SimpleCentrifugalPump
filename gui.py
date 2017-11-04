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
        self.f_tab2.grid(row=0, column=0)

        n.add(self.f_tab1, text="Feasability Study")
        n.add(self.f_tab2, text="Pump Shaft Diameter")


def main():
    root = tk.Tk()

    app = Gui(master=root)
    app.master.title('Dimensioning Radial Pump')
    app.mainloop()

if __name__ == '__main__':
    main()
