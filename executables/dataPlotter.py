import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import FormatStrFormatter

class axial_data_plotter:
    fsize = 15
    ticksize = 13
    P_t_in = 155000.33
    T_t_in = 308.0
    fracfactor = 1.3
    tickformat = '%.2f'
    figsize = (6, 5)
    i_fig = 0

    def __init__(self):
        if not os.path.isdir(os.getcwd()+'/Performance_Data'):
            os.system("mkdir Performance_Data")
        self.direct = os.getcwd() + '/Performance_Data/'
        self.extract_output_data()
        self.plot_data()


    def plot_data(self):
        data_list = list(self.Data.__dict__.keys())
        os.system("mkdir " + self.direct + "Images")
        for data in data_list:
            os.system("mkdir " + self.direct + "Images/" + data)
            if data == 'axial_data':
                data_object = getattr(self.Data, data)
                vars_list = list(data_object.__dict__.keys())
                vars_list.remove('X')
                for var in vars_list:
                    x = data_object.X
                    y = getattr(data_object, var)
                    x_label = 'axial coordinate'
                    y_label = var
                    title = var + ' axial data plot'
                    filename = var
                    self.plot_figure(x=x, y=y, x_label=x_label, y_label=y_label, title=title, directory=self.direct + "Images/" + data, filename=filename)
            else:
                data_object = getattr(self.Data, data)
                vars_list = list(data_object.__dict__.keys())
                vars_list.remove('X')
                for var in vars_list:
                    y = data_object.X
                    x = getattr(data_object, var)
                    x_label = var
                    y_label = 'radial coordinate'
                    title = data + ' ' + var + ' radial data plot'
                    filename = 'radial_'+var
                    self.plot_figure(x=x, y=y, x_label=x_label, y_label=y_label, title=title, directory=self.direct + "Images/" + data, filename=filename)

    def plot_figure(self, x, y, x_label, y_label, title, directory, filename):
        fig = plt.figure(self.i_fig, figsize=self.figsize)
        ax = fig.add_subplot(111)
        plt.title(title, fontsize=self.fsize)
        ax.plot(x, y)
        plt.xlabel(x_label, fontsize=self.fsize)
        plt.ylabel(y_label, fontsize=self.fracfactor * self.fsize)
        ax.grid()
        ax.tick_params(axis='both', which='both', labelsize=self.ticksize)
        ax.xaxis.set_major_formatter(FormatStrFormatter(self.tickformat))
        ax.yaxis.set_major_formatter(FormatStrFormatter(self.tickformat))
        fig.savefig(directory + '/'+filename+".eps", bbox_inches='tight', format='eps')
        plt.close(self.i_fig)
        self.i_fig += 1

    def extract_output_data(self):
        output_files = os.listdir(self.direct)
        output_files.remove("Machine_objectives.txt")
        setattr(self, 'Data', type('', (), {})())
        for file in output_files:
            setattr(self.Data, file[:-4], type('', (), {})())
            station_object = getattr(self.Data, file[:-4])
            flow_data_names = open(self.direct + file, 'r').readline().strip().split(sep='\t')
            flow_data = np.loadtxt(self.direct + file, delimiter='\t', skiprows=1)
            for i in range(len(flow_data_names)):
                setattr(station_object, flow_data_names[i], flow_data[:, i])

