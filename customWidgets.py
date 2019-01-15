import tkinter as tk
import tkinter.ttk as ttk
from PIL import ImageTk, Image


#only pass one of unitvar or unittext, titlevar or titletext
class EntryFrame(): 
    def __init__(self, parent, *args, **kwargs):
        self.frame = tk.Frame(parent)

        #there's definitely a better way to do this
        for key in kwargs:
            if (key=='title'):
                self.title=ttk.Label(self.frame, text=kwargs[key])
            if (key=='unit'):
                self.unit=ttk.Label(self.frame, text=kwargs[key])
            if (key=='titlevar'):
                self.title=ttk.Label(self.frame, textvariable=kwargs[key])
            if (key=='unitvar'):
                self.unit=ttk.Label(self.frame, textvariable=kwargs[key])
            if (key=='entryvar'):
                self.entryVar=kwargs[key]
            if (key=='caption'):
                self.caption=kwargs[key]

        if (kwargs.get('titlevar') is None) and (kwargs.get('title') is None):
            self.title=ttk.Label(self.frame, text='Title')
        if (kwargs.get('unitvar') is None) and (kwargs.get('unit') is None):
            self.unit=ttk.Label(self.frame, text='(unit)')           
        if (kwargs.get('entryvar') is None):
            self.entryVar = tk.StringVar()
        if (kwargs.get('caption') is None):
            self.caption='caption'

        self.caption = ttk.Label(self.frame, text=self.caption)
        self.entry = ttk.Entry(self.frame, textvariable=self.entryVar)
        
        self.title.grid(row=0, column=0, sticky='es')
        self.unit.grid(row=0, column=1, sticky='ws')
        self.entry.grid(row=1, column=0, columnspan=2, padx=3, pady=3, sticky='n')
        self.caption.grid(row=1, column=2, padx=2)
        self.grid=self.frame.grid

class OutputFrame(tk.Frame):
        def __init__(self, parent, *args, **kwargs):
            self.frame = tk.Frame(parent)

            self.Style = ttk.Style()
            self.Style.configure('OF.Tframe', background='#111111')
            self.opContainer = ttk.Frame(self.frame, relief='sunken', height=25, width=135)

                    #there's definitely a better way to do this
            for key in kwargs:
                if (key=='title'):
                    self.title=ttk.Label(self.frame, text=kwargs[key])
                if (key=='unit'):
                    self.unit=ttk.Label(self.opContainer, text=kwargs[key])
                if (key=='titlevar'):
                    self.title=ttk.Label(self.frame, textvariable=kwargs[key])
                if (key=='unitvar'):
                    self.unit=ttk.Label(self.opContainer, textvariable=kwargs[key])
                if (key=='variable'):
                    self.variable=kwargs[key]
                if (key=='caption'):
                    self.caption=kwargs[key]
                

            if (kwargs.get('titlevar') is None) and (kwargs.get('title') is None):
                self.title=ttk.Label(self.frame, text='Title')
            if (kwargs.get('unitvar') is None) and (kwargs.get('unit') is None):
                self.unit=ttk.Label(self.opContainer, text='(unit)')           
            if (kwargs.get('caption') is None):
                self.caption=''
            if (kwargs.get('variable') is None):
                self.variable=tk.StringVar()
                self.variable.set('')
                #print(self.variable.get())

            self.caption = ttk.Label(self.frame, text=self.caption)
            
            #print(self.variable.get())
            #WHY THE FUCK WON'T YOU OUTPUT
            self.outputvar = tk.Label(self.opContainer, textvariable=self.variable)


            self.title.grid(row=0, column=0, sticky='es')
            self.unit.grid(row=0, column=1, padx=3, pady=3, sticky='nsw')
            self.opContainer.grid(row=1, column=0, columnspan=2, padx=3, pady=3, sticky='n')
            self.opContainer.grid_propagate(0)
            self.caption.grid(row=1, column=2, padx=2)
            self.outputvar.grid(row=0, column=0, padx=3, pady=3, sticky='nsw')
                                         
            self.grid=self.frame.grid
        

class TripleEntryFrame(): 
    def __init__(self, parent, *args, **kwargs):
        self.frame = tk.Frame(parent)

        #there's definitely a better way to do this
        for key in kwargs:
            if (key=='title'):
                self.title=ttk.Label(self.frame, text=kwargs[key])
            if (key=='unit'):
                self.unit=ttk.Label(self.frame, text=kwargs[key])
            if (key=='titlevar'):
                self.title=ttk.Label(self.frame, textvariable=kwargs[key])
            if (key=='unitvar'):
                self.unit=ttk.Label(self.frame, textvariable=kwargs[key])
            if (key=='caption'):
                self.caption=kwargs[key]
            if (key=='xvariable'):
                self.xentry = ttk.Entry(self.frame, textvariable=kwargs[key])
            if (key=='yvariable'):
                self.yentry = ttk.Entry(self.frame, textvariable=kwargs[key])
            if (key=='zvariable'):
                self.zentry = ttk.Entry(self.frame, textvariable=kwargs[key])

        if (kwargs.get('titlevar') is None) and (kwargs.get('title') is None):
            self.title=ttk.Label(self.frame, text='Title')
        if (kwargs.get('unitvar') is None) and (kwargs.get('unit') is None):
            self.unit=ttk.Label(self.frame, text='(unit)')           
        if (kwargs.get('caption') is None):
            self.caption=''
        if (kwargs.get('xvariable') is None):
            self.xentry = ttk.Entry(self.frame)
        if (kwargs.get('yvariable') is None):
            self.yentry = ttk.Entry(self.frame)
        if (kwargs.get('zvariable') is None):
            self.zentry = ttk.Entry(self.frame)
        

        self.caption = ttk.Label(self.frame, text=self.caption)
    
        self.xlabel = ttk.Label(self.frame, text='X')
        self.ylabel = ttk.Label(self.frame, text='Y')
        self.zlabel = ttk.Label(self.frame, text='Z')
        
        self.title.grid(row=0, column=1, sticky='es')
        self.unit.grid(row=0, column=2, sticky='ws')
        self.xentry.grid(row=1, column=1, columnspan=2, padx=3, pady=3, sticky='n')
        self.yentry.grid(row=2, column=1, columnspan=2, padx=3, pady=3, sticky='n')
        self.zentry.grid(row=3, column=1, columnspan=2, padx=3, pady=3, sticky='n')
        self.xlabel.grid(row=1, column=0, padx=3, pady=3, sticky='e')
        self.ylabel.grid(row=2, column=0, padx=3, pady=3, sticky='e')
        self.zlabel.grid(row=3, column=0, padx=3, pady=3, sticky='e')
        self.grid=self.frame.grid
        
        
class CalculationFrame(tk.Frame):
    def __init__ (self, parent, *args, **kwargs):
        tk.Frame.__init__(self, parent=None, *args, **kwargs)
        self.parent = parent


        self.imgContainer = ttk.Frame(parent, relief='groove')
        self.imgContainer.grid(row=0, column=0, padx=10, pady=10)
 
        self.imgFrame = ttk.Frame(self.imgContainer, relief='groove')
        self.imgFrame.grid(row=0, column=0, padx=10, pady=10)
    
        self.captionFrame = ttk.Frame(self.imgContainer)
        self.captionFrame.grid(row=1, column=0, padx=10, pady=10, sticky='nwe')



        self.calcsContainer = ttk.Frame(parent, relief='groove')
        self.calcsContainer.grid(row=0, column=1, padx=10, pady=10)

        self.inputsFrame = ttk.Frame(self.calcsContainer, relief='groove')
        self.inputsFrame.grid(row=0, column=0, padx=10, pady=10, ipadx=7, ipady=5)
        ipTitle = ttk.Label(self.inputsFrame, font=12, text='Inputs')
        ipTitle.grid(row=0, column=0)
        
        self.outputsFrame = ttk.Frame(self.calcsContainer, relief='groove')
        self.outputsFrame.grid(row=1, column=0, padx=10, pady=10, ipadx=7, ipady=5)
        opTitle = ttk.Label(self.outputsFrame, font=12, text='Outputs')
        opTitle.grid(row=0, column=0)
    

class TestApplication(tk.Frame):
    def __init__(self, master=None, *args, **kwargs):
        tk.Frame.__init__(self, master=None, *args, **kwargs)
        self.master=master
        self.grid(row=0, column=0)

        calcFrameTest = CalculationFrame(self)
        entryFrameTest = TripleEntryFrame(self)
        entryFrameTest.grid(row=1,column=1)


        
        

        
        
        


if __name__ == "__main__":
    root = tk.Tk()
    Applicatiion = TestApplication(root)
    root.mainloop()
