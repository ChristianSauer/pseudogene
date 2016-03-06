from Tkconstants import END
import os
import tkMessageBox
import ttk
from pseudogene import PseudoGeneFinder

__author__ = 'Christian'
import Tkinter
from tkFileDialog import askopenfilename


class PseudoGeneGUI(Tkinter.Tk):
    def __init__(self, parent):
        Tkinter.Tk.__init__(self, parent)
        self.initialize()
        self.parent = parent
        self.filename = None

    def initialize(self):
        self.grid()

        row = 0
        label_txt = "Choose Config:"
        label = Tkinter.Label(self, anchor="w", text=label_txt)
        label.grid(column=0, row=0, columnspan=2, sticky='EW')

        button = Tkinter.Button(self, text=u"Choose a config.json", command=self.on_file_chooser, anchor="e")
        button.grid(column=1, row=0)

        row += 1
        label_txt = "Genes to search for. Separate genes into multiple Lines"
        label = Tkinter.Label(self, anchor="w", text=label_txt)
        label.grid(column=0, row=1, columnspan=2, sticky='EW')

        row += 1
        # noinspection PyAttributeOutsideInit In itialize called before init?!?
        self.entry = Tkinter.Text(self)
        self.entry.grid(column=0, row=2, columnspan=2, rowspan=20, sticky='EW')

        row += 1
        label_txt = "Execute Pseudogene search"
        label = Tkinter.Label(self, anchor="w", text=label_txt)
        label.grid(column=0, row=24, columnspan=2, sticky='EW')

        button = Tkinter.Button(self, text=u"Execute", command=self.on_execute, anchor="e")
        button.grid(column=1, row=24)

        # self.progress = ttk.Progressbar(self, orient="horizontal", length=200, mode="indeterminate")
        # self.progress.grid(column=0, row=25, columnspan=3)

    def on_file_chooser(self):
        filename = askopenfilename(defaultextension="*.json", initialdir=os.getcwd(), multiple=False,
                                   parent=self.parent,
                                   title="Please choose a config.json", filetypes=[('json files', '.json')])

        if filename is None:
            tkMessageBox.showerror("Error", "Please choose a valid config.json")

        self.filename = filename

    def on_execute(self):
        genes = self.entry.get("1.0", END)
        if self.filename is None or genes is None:
            tkMessageBox.showerror("Error", "Please provide a valid config.json and some genes")

        finder = PseudoGeneFinder(self.filename, genes)
        # self.progress.start()
        finder.find_pseudogenes()
        tkMessageBox.showinfo("OK", "Found all requested pseudogenes")
        # self.progress.stop()


if __name__ == "__main__":
    app = PseudoGeneGUI(None)
    app.title('PseudoGene Finder')
    app.mainloop()
