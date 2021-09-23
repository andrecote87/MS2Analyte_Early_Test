import sys
import os
from PyQt5 import QtCore, QtWidgets
from qtmodern.styles import dark
from qtmodern.windows import ModernWindow
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
import pyqtgraph as pg
import numpy as np
import pandas
import random
import time
import pickle
from pyqtgraph.Qt import QtCore, QtGui
from decimal import Decimal, ROUND_DOWN
import cProfile
import pstats
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QLabel, QLineEdit, QTabWidget, \
				  QGridLayout, QVBoxLayout, QHBoxLayout, QGroupBox, QDialog
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot

import pyqtgraph.console
from collections import namedtuple
from itertools import chain
import glob
from pathlib import Path
from ms2analyte.file_handling import data_import
from ms2analyte.file_handling import file_load
from ms2analyte.visualizations import file_open_model
from pyteomics import mgf, auxiliary
import numpy
from matchms import Spectrum
from matchms.exporting import save_as_mgf
from zipfile import ZipFile
import ms2analyte.file_open_dialogue
from threading import *
from PyQt5.QtWidgets import *
from ms2analyte.file_open_dialogue import MainWindowUIClass
from PyQt5 import uic
from PyQt5.QtWidgets import QMessageBox


global open_state
open_state = 1
# Handle high resolution displays:
if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)




########### Data Import #################
lastClicked = []
clickedPen = pg.mkPen('s', width=5)

clickedPen_legend = pg.mkPen('s', width=3)
lastClicked_legend = []





def fetch_sample_list(input_structure, input_type):

	if open_state ==2:
		sample_list = data_import.name_extract(input_structure, input_type)
	else:
		sample_list = []

	return sample_list



#######################################################################################
def import_dataframe(input_structure, input_type, sample_name):
	print('TEST importdataframe 1')
	if open_state == 2:
		with open((os.path.join(output_path, experiment_name + "_experiment_import_parameters.pickle")), 'rb') as f:
			experiment_info = pickle.load(f)



		print(experiment_info.blanks_exist)
		print('TEST importdataframe 2')

		print(input_structure)
		print(sample_name)
		print(input_type)
		if experiment_info.blanks_exist == True:
			with open((os.path.join(input_structure, input_type, sample_name + "_all_replicates_blanked_dataframe.pickle")), 'rb') as f:
				print('TEST importdataframe 3')
				df = pickle.load(f)
				# df = pandas.read_pickle(f)

		else:

			with open((os.path.join(input_structure, input_type, sample_name + "_all_replicates_dataframe.pickle")), 'rb') as f:
				print('TEST importdataframe 3')
				df = pickle.load(f)
			print('TEST importdataframe 4')

		print('TEST importdataframe 5')
		print(df)
	else:
		df = []
	print('TEST importdataframe 6')
	return df



#######################################################################################
def import_experiment_dataframe(input_structure):
	if open_state == 2:
		with open((os.path.join(input_structure,  experiment_name + "_experiment_analyte_overview_tableau_output.csv"))) as f:
			df = pandas.read_csv(f)
	else:
		df = []
	return df
    




#######################################################################################
def import_ms1_dataframe(input_structure, input_type, sample_name):


	if open_state == 2:
		print('test1')
		with open((os.path.join(input_structure, input_type, sample_name + '_replicate_analyte_mass_spectra.pickle')), 'rb') as f:


			inputdata = pickle.load(f)
			
		print('test2')
		print('ms1 inputdata',inputdata)
		b = 0
		replicate_analyte_id = []
		average_mass = []
		relative_intensity = []

		while b < len(inputdata):
			i = 0
			while i <len(inputdata[b].replicate_analyte_ms1_spectrum):
				replicate_analyte_id.append(inputdata[b].replicate_analyte_id)
				average_mass.append(inputdata[b].replicate_analyte_ms1_spectrum[i].average_mass)
				relative_intensity.append(inputdata[b].replicate_analyte_ms1_spectrum[i].relative_intensity)
				i+=1


			b+=1
		np.round(average_mass,1)
		relative_intensity_zeros = [0]*len(relative_intensity)
		average_mass_lower = [0]*len(average_mass)
		i=0
		while i < len(average_mass):
			average_mass_lower[i] = average_mass[i] - 0.0001
			i+=1

		average_mass_upper = [0]*len(average_mass)
		i=0
		while i < len(average_mass):
			average_mass_upper[i] = average_mass[i] + 0.0001
			i+=1

		data1 = {'replicate_analyte_id': replicate_analyte_id,
				'average_mass': average_mass,
				'relative_intensity': relative_intensity
				}
		data2 = {'replicate_analyte_id': replicate_analyte_id,
				'average_mass': average_mass_upper,
				'relative_intensity': relative_intensity_zeros
				}

		data3 = {'replicate_analyte_id': replicate_analyte_id,
				'average_mass': average_mass_lower,
				'relative_intensity': relative_intensity_zeros
				}

		df1 = pandas.DataFrame (data1, columns = ['replicate_analyte_id', 'average_mass','relative_intensity'])

		df2 = pandas.DataFrame (data2, columns = ['replicate_analyte_id', 'average_mass','relative_intensity'])

		df3 = pandas.DataFrame (data3, columns = ['replicate_analyte_id', 'average_mass','relative_intensity'])


		df_combine = [df1,df2,df3]

		df_combine = pandas.concat(df_combine)



		print('test3')
		with open((os.path.join(input_structure, input_type, sample_name + '_R1_analytes.pickle')), 'rb') as f:
			inputdata2 = pickle.load(f)
		print('test4')
		b = 0
		analyte_id = []
		max_peak_intensity_mass = []
		max_peak_intensity = []


		while b < len(inputdata2):

			analyte_id.append(inputdata2[b].analyte_id)
			max_peak_intensity_mass.append(inputdata2[b].max_peak_intensity_mass)
			max_peak_intensity.append(inputdata2[b].max_peak_intensity)



			b+=1


		max_peak_intensity_zeros = [0]*len(max_peak_intensity)
		max_peak_intensity_mass_lower = [0]*len(max_peak_intensity_mass)
		i=0
		while i < len(max_peak_intensity_mass):
			max_peak_intensity_mass_lower[i] = max_peak_intensity_mass[i] - 0.0001
			i+=1





		max_peak_intensity_zeros = [0]*len(max_peak_intensity)
		max_peak_intensity_mass_upper = [0]*len(max_peak_intensity_mass)
		i=0
		while i < len(max_peak_intensity_mass):
			max_peak_intensity_mass_upper[i] = float(max_peak_intensity_mass[i] + 0.0001)
			i+=1


		data4 = {'analyte_id': analyte_id,
				'max_peak_intensity_mass': max_peak_intensity_mass,
				'max_peak_intensity': max_peak_intensity
				}
		data5 = {'analyte_id': analyte_id,
				'max_peak_intensity_mass': max_peak_intensity_mass_upper,
				'relative_intensity': max_peak_intensity_zeros
				}

		data6 = {'analyte_id': analyte_id,
				'max_peak_intensity_mass': max_peak_intensity_mass_lower,
				'max_peak_intensity': max_peak_intensity_zeros
				}

		df4 = pandas.DataFrame (data4, columns = ['analyte_id', 'max_peak_intensity_mass','max_peak_intensity'])

		df5 = pandas.DataFrame (data5, columns = ['analyte_id', 'max_peak_intensity_mass','max_peak_intensity'])

		df6 = pandas.DataFrame (data6, columns = ['analyte_id', 'max_peak_intensity_mass','max_peak_intensity'])

		df_combine2 = [df4,df5,df6]

		df_combine2 = pandas.concat(df_combine2)







		print('test5')
		with open((os.path.join(input_structure,input_type,  experiment_name + '_experiment_analyte_mass_spectra.pickle')), 'rb') as f:




			inputdata3 = pickle.load(f)
			
		print('test6')
		b = 0
		experiment_analyte_id = []
		average_mass = []
		relative_intensity = []

		while b < len(inputdata3):
			i = 0
			while i <len(inputdata3[b].relative_experiment_mass_spectrum):
				experiment_analyte_id.append(inputdata3[b].experiment_analyte_id)
				average_mass.append(inputdata3[b].relative_experiment_mass_spectrum[i].average_mass)
				relative_intensity.append(inputdata3[b].relative_experiment_mass_spectrum[i].relative_intensity)
				i+=1


			b+=1
		np.round(average_mass,1)
		relative_intensity_zeros = [0]*len(relative_intensity)
		average_mass_lower = [0]*len(average_mass)
		i=0
		while i < len(average_mass):
			average_mass_lower[i] = average_mass[i] - 0.0001
			i+=1

		average_mass_upper = [0]*len(average_mass)
		i=0
		while i < len(average_mass):
			average_mass_upper[i] = average_mass[i] + 0.0001
			i+=1

		data7 = {'experiment_analyte_id': experiment_analyte_id,
				'average_mass': average_mass,
				'relative_intensity': relative_intensity
				}
		data8 = {'experiment_analyte_id': experiment_analyte_id,
				'average_mass': average_mass_upper,
				'relative_intensity': relative_intensity_zeros
				}

		data9 = {'experiment_analyte_id': experiment_analyte_id,
				'average_mass': average_mass_lower,
				'relative_intensity': relative_intensity_zeros
				}

		df7 = pandas.DataFrame (data7, columns = ['experiment_analyte_id', 'average_mass','relative_intensity'])

		df8 = pandas.DataFrame (data8, columns = ['experiment_analyte_id', 'average_mass','relative_intensity'])

		df9 = pandas.DataFrame (data9, columns = ['experiment_analyte_id', 'average_mass','relative_intensity'])


		df_combine3 = [df7,df8,df9]

		df_combine3 = pandas.concat(df_combine3)
		print('test7')
	else:
		df_combine = []
		df_combine2 = []
		df_combine3 = []

	return df_combine, df_combine2, df_combine3



def import_ms2_dataframe(input_structure, input_type, sample_name):

	if open_state == 2:

		with open((os.path.join(output_path, experiment_name + "_experiment_import_parameters.pickle")), 'rb') as f:
			experiment_info = pickle.load(f)






		if experiment_info.ms2_type == 'DIA':

			with open((os.path.join(output_path, input_type, sample_name + '_replicate_analyte_mass_spectra.pickle')), 'rb') as f:


				inputdata = pickle.load(f)
				
			print('ms2 inputdata ****',inputdata)
			b = 0
			replicate_analyte_id = []
			average_mass = []
			relative_intensity = []


			
			while b < len(inputdata):




				if inputdata[b].replicate_analyte_ms2_spectrum is not None:
					length_ms2 = len(inputdata[b].replicate_analyte_ms2_spectrum)
					i = 0
					while i < length_ms2:


						replicate_analyte_id.append(inputdata[b].replicate_analyte_id)
						average_mass.append(inputdata[b].replicate_analyte_ms2_spectrum[i].average_mass)
						relative_intensity.append(inputdata[b].replicate_analyte_ms2_spectrum[i].relative_intensity)
						i+=1
				else:
					b+=1


				b+=1

			print('REPLICATE ANALYTE ID MS2',replicate_analyte_id)
			print('AVERAGE MASS MS2',average_mass)
			print('RELATIVE INTENSITY MS2',relative_intensity)

			np.round(average_mass,1)
			relative_intensity_zeros = [0]*len(relative_intensity)
			average_mass_lower = [0]*len(average_mass)
			i=0
			while i < len(average_mass):
				average_mass_lower[i] = average_mass[i] - 0.0001
				i+=1

			average_mass_upper = [0]*len(average_mass)
			i=0
			while i < len(average_mass):
				average_mass_upper[i] = average_mass[i] + 0.0001
				i+=1

			data1 = {'replicate_analyte_id': replicate_analyte_id,
					'average_mass': average_mass,
					'relative_intensity': relative_intensity
					}
			data2 = {'replicate_analyte_id': replicate_analyte_id,
					'average_mass': average_mass_upper,
					'relative_intensity': relative_intensity_zeros
					}

			data3 = {'replicate_analyte_id': replicate_analyte_id,
					'average_mass': average_mass_lower,
					'relative_intensity': relative_intensity_zeros
					}

			df1 = pandas.DataFrame (data1, columns = ['replicate_analyte_id', 'average_mass','relative_intensity'])

			df2 = pandas.DataFrame (data2, columns = ['replicate_analyte_id', 'average_mass','relative_intensity'])

			df3 = pandas.DataFrame (data3, columns = ['replicate_analyte_id', 'average_mass','relative_intensity'])


			df_combine = [df1,df2,df3]

			df_combine = pandas.concat(df_combine)

			print('MS2 DF COMBINE',df_combine)
		if experiment_info.ms2_type == 'DDA':




			with open((os.path.join(output_path,input_type,  sample_name + '_R1_analytes.pickle')), 'rb') as f:


				inputdata = pickle.load(f)




			Total_Analytes = len(inputdata)



			Analyte_List = []

			ms1_List = []

			ms1_rt_List = []

			ms2_List = []


			Intensity_List = []


			i = 0
			while i <  Total_Analytes:

				b = 0
				Total_Av_Mass = len(inputdata[i].peak_list)
				while b < Total_Av_Mass:

					if inputdata[i].peak_list[b].dda_data is not None:
						c = 0
						DDA_Length = len(inputdata[i].peak_list[b].dda_data[0])

					
						while c < DDA_Length:

							Analyte_List.append(inputdata[i].analyte_id)
							ms1_List.append(inputdata[i].peak_list[b].average_mass)
							ms1_rt_List.append(inputdata[i].peak_list[b].rt)
							ms2_List.append(inputdata[i].peak_list[b].dda_data[0][c])
							Intensity_List.append(inputdata[i].peak_list[b].dda_data[1][c])
							c+=1
					else:
						pass
					b+=1
				i+=1





			intensity_zeros = [0]*len(Intensity_List)





			ms2_lower = [0]*len(ms2_List)
			i=0
			while i < len(ms2_List):
				ms2_lower[i] = ms2_List[i] - 0.000001
				i+=1

			ms2_upper = [0]*len(ms2_List)
			i=0
			while i < len(ms2_List):
				ms2_upper[i] = ms2_List[i] + 0.000001
				i+=1


			data1 = {'analyte_id': Analyte_List,
					'ms1_average_mass': ms1_List,
					'ms1_rt': ms1_rt_List,
					'ms2_data': ms2_List,
					'Intensity': Intensity_List
					}


			data2 = {'analyte_id': Analyte_List,
					'ms1_average_mass': ms1_List,
					'ms1_rt': ms1_rt_List,
					'ms2_data': ms2_lower,
					'Intensity': intensity_zeros
					}


			data3 = {'analyte_id': Analyte_List,
					'ms1_average_mass': ms1_List,
					'ms1_rt': ms1_rt_List,
					'ms2_data': ms2_upper,
					'Intensity': intensity_zeros
					}

			df1 = pandas.DataFrame (data1, columns = ['analyte_id', 'ms1_average_mass','ms1_rt','ms2_data','Intensity'])
			df2 = pandas.DataFrame (data2, columns = ['analyte_id', 'ms1_average_mass','ms1_rt','ms2_data','Intensity'])
			df3 = pandas.DataFrame (data3, columns = ['analyte_id', 'ms1_average_mass','ms1_rt','ms2_data','Intensity'])

			print('*******************************************************************',df1[df1.analyte_id == 7])
			df_combine = [df1,df2,df3]

			df_combine = pandas.concat(df_combine)


			print('*******************************************************************',df_combine[df_combine.analyte_id == 7])


	else:
		df_combine = []





	return df_combine
	
#################################################################################################
#### Fetch Input Structure, Input type and Sample name list

# input_structure = data_import.input_data_structure()


class NotEmptyValidator(QValidator):
	def validate(self, text:str, pos):
		if bool(text.strip()):
			state = QValidator.Acceptable
		else:
			state = QValidator.Intermediate
		return state, text, pos



def getFiles(path):
	answer = []
	endPaths = glob.glob(path + "*parameters.pickle")
	answer += endPaths
	if len(glob.glob(path + "*Output/")) > 0:
		answer += getFiles(path + "*Output/")
	return answer
filepaths = getFiles("./")
print('********************************************** ANSWER',filepaths)







# input_structure = file_load.sample_dataframe_concat()

# print(input_structure)















input_type = "Samples"


# path_name_list = []

# for path in Path(data_import.input_data_structure().output_directory).rglob('*parameters.pickle'):
	
# 	path_name_list.append(path.name)
# path_name = str(path_name_list[0])
# print('***********************************************************************************PATH NAME',path_name_list)




# object1 = pandas.read_pickle((os.path.join(data_import.input_data_structure().output_directory, path_name)))

# input_structure = object1

# print('************************************************ experiment name',object1.experiment_name)


# experiment_name = object1.experiment_name




# print('***********************************************OUTPUT DIRECTORY',object1.output_directory)


# objects = []
# with (open(os.path.join(input_structure.output_directory,  'ms2_experiment_experiment_import_parameters.pickle'), "rb")) as openfile:
#     while True:
#         try:
#             objects.append(pickle.load(openfile))
#         except EOFError:
#             break
# with open((os.path.join(input_structure.output_directory, input_type, 'ms2_experiment_experiment_import_parameters.pickle'), 'rb')) as f:

# 	input_data = pickle.load(f)




# sample_names = fetch_sample_list(object1, input_type)

# print('SAMPLE NAMES',sample_names)


#################################### Start GUI  ##################################################


class State:
	def __init__(self, null, blank, massMin, massMax, rtMin, rtMax, Toggle_rt):

		self.null = null
		self.blank = blank
		self.massMin = massMin
		self.massMax = massMax
		self.rtMin = rtMin
		self.rtMax = rtMax
		self.Toggle_rt = Toggle_rt


	def return_null_state(self):


		return self.null


	def return_blank_state(self):


		return self.blank


	def return_massMin_state(self):

		
		return self.massMin



	def return_massMax_state(self):

		
		return self.massMax


	def return_rtMin_state(self):

		# print("State:", self.rtMin)
		return self.rtMin


	def return_rtMax_state(self):

		# print("State:", self.rtMax)
		return self.rtMax

	def return_Toggle_rt(self):

		return self.Toggle_rt


class Tab_State:
	def __init__(self, sample_state, replicate_state, experiment_state, diversity_state):

		self.sample_state = sample_state
		self.replicate_state = replicate_state
		self.experiment_state = experiment_state
		self.diversity_state = diversity_state

	def return_sample_state(self):


		return self.sample_state

	def return_replicate_state(self):


		return self.replicate_state

	def return_experiment_state(self):

		# print("Experiment State:", self.experiment_state)
		return self.experiment_state

	def return_diversity_state(self):

		# print("Diversity State:", self.diversity_state)

		return self.diversity_state








class Current_Sample:
	def __init__(self, current_sample):


		self.current_sample = current_sample


	def return_current_sample(self):

		print("CURRENT SAMPLE", self.current_sample)

		return self.current_sample












class Window (QDialog):



	def __init__(self):
		super(Window, self).__init__()
		self.right_button_is_clicked = 0


        # changing the background color to yellow 

		
		app = QtWidgets.QApplication(sys.argv)
		#dark(app) 
		app.setStyle('Fusion')
		# self.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")

		self.setWindowTitle("MS2Analyte")
		self.setGeometry(100,100,1700,900)
		self.setWindowIcon(QIcon('ms2analyte_icon.png'))
		pg.setConfigOption('background', 'w')
		pg.setConfigOption('foreground', 'k')


		self.region_1()
		self.Plot_Sample()
		self.Plot_ms1()
		self.Plot_ms2()
		self.Plot_Analyte_Legend()
		windowLayout = QVBoxLayout()
		# windowLayout.setStyleSheet("background-color: white;")


		self.menu_bar = QMenuBar()





		self.file_menu = self.menu_bar.addMenu('File')
		self.analysis_menu = self.menu_bar.addMenu('Analysis')
		self.preferences_menu = self.menu_bar.addMenu('Preferences' )
		self.export_menu = self.menu_bar.addMenu('Export')
		self.help_menu = self.menu_bar.addMenu('Help')





		self.new_analysis_action = QAction('New Analysis',app)
		self.analysis_menu.addAction(self.new_analysis_action)
		self.new_action = QAction('New', app)
		self.open_action = QAction('Open', app)
		self.open_recent_action = QAction('Open Recent', app)
		self.exit_action = QAction('Exit', app)
		self.exit_action.triggered.connect(exit)

		self.file_menu.addAction(self.open_action)

		self.file_menu.addAction(self.exit_action)



		self.readme_action = QAction('Read Me', app)






		self.help_menu.addAction(self.readme_action)
		self.readme_action.triggered.connect(self.read_me)

		self.open_action.triggered.connect(self.openFileNamesDialog)
		self.open_action.triggered.connect(self.info_test)

		self.open_action.triggered.connect(self.Plot_Sample)

		self.open_action.triggered.connect(self.Plot_Analyte_Legend)



		self.new_analysis_action.triggered.connect(self.launch_analysis)


		# self.open_action.triggered.connect(self.region_1)


		# self.open_action.triggered.connect(self.Plot_Sample)



		self.csv_action = QAction('Export CSV', app)
		self.export_menu.addAction(self.csv_action)
		self.csv_action.triggered.connect(self.sample_export)

		self.mgf_action = QAction('Export MGF', app)
		self.export_menu.addAction(self.mgf_action)
		self.mgf_action.triggered.connect(self.ms1_DF_Export)

		windowLayout.addWidget(self.menu_bar)
		windowLayout.addWidget(self.horizontalGroupBox)
		self.setLayout(windowLayout)

 

		# self.setStyleSheet("background-color: white;")


		
######## Hide some plots initially ########
		
		self.mz_r1_plot.hide()
		self.mz_r2_plot.hide()
		self.mz_r3_plot.hide()

		self.rt_r1_plot.hide()
		self.rt_r2_plot.hide()
		self.rt_r3_plot.hide()
		

		self.ex_bOn_plot.hide()
		self.div_plot.hide()


		# self.showMaximized()
		self.setWindowFlags(self.windowFlags() |
			QtCore.Qt.WindowMinimizeButtonHint |
			QtCore.Qt.WindowSystemMenuHint)

		self.setWindowFlags(self.windowFlags() |
			QtCore.Qt.WindowMaximizeButtonHint |
			QtCore.Qt.WindowSystemMenuHint)



		self.show()
		start_message = QMessageBox()
		start_message.setWindowTitle("MS2Analyte")
		start_message.setText('Welcome to MS2Analyte. To get started, run a new analysis or open an existing output folder')
		x = start_message.exec_()
	def sample_export(self):
		state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), self.rtMinBox.text(), self.rtMaxBox.text(),self.Toggle_rt.isChecked())



		BlankState = state.return_blank_state()
		NullState = state.return_null_state()
		comboText = self.comboBox.currentText()


		sorted_df_mz, sorted_df_rt= self.Sample_Plot_DF()





		if BlankState or NullState == True:
			name, _ = QFileDialog.getSaveFileName(self, "Save", os.getcwd()+os.sep+comboText +'_Filter', "CSV Files (*.csv)")

		else:
			name, _ = QFileDialog.getSaveFileName(self, "Save", os.getcwd()+os.sep+comboText, "CSV Files (*.csv)")

		df = sorted_df_mz.to_csv()

		if not name : return 0

		sorted_df_mz.to_csv(name, index=False)




	def region_1(self):
		self.horizontalGroupBox = QGroupBox()




		layoutMain = QGridLayout()


		# self.fileMenuBar = QMenuBar(self).addMenu('File')



		# self.fileMenuBar.addMenu('File')
		# self.preferencesMenuBar = QMenuBar(self)
		# self.preferencesMenuBar.addMenu('Preferences')

		# self.s1 = QScrollBar()

		# self.s1.sliderMoved.connect(self.Plot_Analyte_Legend)


		self.btn2 = QPushButton("Sample", self)

		self.btn2.clicked.connect(self.loading_sample)
		self.btn2.clicked.connect(self.Sample_Plot_DF)
		self.btn2.clicked.connect(self.show_sample)
		self.btn2.clicked.connect(self.Plot_Sample)
		self.btn2.clicked.connect(self.hide_replicate)
		self.btn2.clicked.connect(self.hide_experiment)
		self.btn2.clicked.connect(self.hide_diversity)

		# self.btn2.resize(buttonX,buttonY)
		# self.btn2.move(20 + buttonmoveX,100 + buttonmoveY)
		self.btn2.setDisabled(True)
		self.btn2.setStyleSheet('background-color : lightgreen')

		
		self.btn3 = QPushButton("Replicate", self)
		self.btn3.clicked.connect(self.loading_replicate)
		self.btn3.clicked.connect(self.Replicate_Plot_DF)
		self.btn3.clicked.connect(self.show_replicate)
		self.btn3.clicked.connect(self.Plot_Replicate)
		self.btn3.clicked.connect(self.hide_sample)

		self.btn3.clicked.connect(self.hide_experiment)
		self.btn3.clicked.connect(self.hide_diversity)


		# self.btn3.resize(buttonX,buttonY)
		# self.btn3.move(170 + buttonmoveX,100 + buttonmoveY)

		self.btn4 = QPushButton("Experiment", self)
		self.btn4.clicked.connect(self.loading_experiment)
		self.btn4.clicked.connect(self.Diversity_Plot_DF)
		self.btn4.clicked.connect(self.show_experiment)
		self.btn4.clicked.connect(self.Plot_Experiment)
		self.btn4.clicked.connect(self.hide_sample)
		self.btn4.clicked.connect(self.hide_replicate)
		self.btn4.clicked.connect(self.hide_diversity)
		# self.btn4.resize(buttonX,buttonY)
		# self.btn4.move(320 + buttonmoveX,100 + buttonmoveY)
		
		self.btn5 = QPushButton("Diversity", self)
		self.btn5.clicked.connect(self.loading_diversity)
		self.btn5.clicked.connect(self.Diversity_Plot_DF)
		self.btn5.clicked.connect(self.show_diversity)
		self.btn5.clicked.connect(self.Plot_Diversity)
		self.btn5.clicked.connect(self.hide_sample)
		self.btn5.clicked.connect(self.hide_replicate)
		self.btn5.clicked.connect(self.hide_experiment)



		# self.btn5.resize(buttonX,buttonY)
		# self.btn5.move(470 + buttonmoveX,100 + buttonmoveY)
		

		
		extractAction = QAction(QIcon('zoom.png'), 'Quit', self)
		extractAction.triggered.connect(self.close_application)
		


########## Grid Region 1 #######################

		layoutRegion1 = QGridLayout()
		layoutRegion1.addWidget(self.btn2,0,1)

		layoutRegion1.addWidget(self.btn3,0,2)
		layoutRegion1.addWidget(self.btn4,0,3)
		layoutRegion1.addWidget(self.btn5,0,4)




##################### Right hand side of interface (Region 2) #########################
		self.comboBox = QComboBox(self)
		self.comboAnalyte = QComboBox(self)
		self.comboBox.currentTextChanged.connect(self.reset_all)
		self.comboBox.currentTextChanged.connect(self.Sample_Plot_DF)
		self.comboBox.currentTextChanged.connect(self.Plot_Sample)
		self.comboBox.currentTextChanged.connect(self.Replicate_Plot_DF)
		self.comboBox.currentTextChanged.connect(self.Plot_Replicate)
		self.comboBox.currentTextChanged.connect(self.Plot_ms1)
		self.comboBox.currentTextChanged.connect(self.Plot_ms2)



		self.chooseSample = QLabel("Choose Sample", self)

		
		# self.chooseAnalyte = QLabel("Choose Analyte", self)

		self.chooseAnalyte = QComboBox(self)
		# self.chooseAnalyte.currentTextChanged.connect(self.fill_combobox)
		self.chooseAnalyte.currentTextChanged.connect(self.Sample_Plot_DF)
		self.chooseAnalyte.currentTextChanged.connect(self.Plot_Sample)

		self.chooseAnalyte.currentTextChanged.connect(self.Replicate_Plot_DF)
		self.chooseAnalyte.currentTextChanged.connect(self.Plot_Replicate)
		self.chooseAnalyte.currentTextChanged.connect(self.ms1_Plot_DF)
		self.chooseAnalyte.currentTextChanged.connect(self.Plot_ms1)
		self.chooseAnalyte.currentTextChanged.connect(self.ms2_Plot_DF)
		self.chooseAnalyte.currentTextChanged.connect(self.Plot_ms2)
		self.chooseAnalyte.currentTextChanged.connect(self.Plot_Analyte_Legend)
		self.chooseAnalyte.addItem('Sample Analyte ID')
		self.chooseAnalyte.addItem('Replicate Analyte ID')
		self.chooseAnalyte.addItem('Experiment Analyte ID')
		self.showBlanks = QLabel("Show Blank Data", self)

		
		self.showNull = QLabel("Show Null Data", self)
		tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
		Sample_State = tab_state.return_sample_state()
		Replicate_State = tab_state.return_replicate_state()




		# if open_state ==2:

		# 	print('======================================= OPEN STATE',open_state)
		# 	print('SAMPLE NAMES******************************************',sample_names)
		# 	n = 0
		# 	while n < len(sample_names):
		# 		self.comboBox.addItem(sample_names[n])
				
		# 		n+=1	




		comboText=self.comboBox.currentText()
		analyteID=self.chooseAnalyte.currentText()
		if open_state ==2:
			print('region 1 TEST TEST ')
			total_df = import_dataframe(output_path, input_type, comboText)
			print('******************************************************************************************* TOTAL DF',total_df)
			# sorted_df = total_df['analyte_id'].replace('',np.nan, inplace=True)
			# sorted_df = total_df.dropna(subset=['analyte_id'],inplace=True)
			# print(sorted_df)


			replicate_id_df = total_df.dropna(subset=['replicate_analyte_id'])
			sorted_df = total_df.dropna(subset=['analyte_id'])
			sorted_df=total_df[total_df.analyte_id != None]
			sorted_df = sorted_df.sort_values(by=['analyte_id'])
			# total_df['analyte_id'].replace('',np.nan, inplace=True)
			mz_array=[]
			analyte_array=[]

			for analyte_id in sorted_df.analyte_id.unique():
				mz_array.append([sorted_df[sorted_df["analyte_id"] == analyte_id]["mz"].to_numpy(),
			                     sorted_df[sorted_df["analyte_id"] == analyte_id]["intensity"].to_numpy()])

				analyte_array.append(analyte_id)


			replicate_id_df = replicate_id_df.dropna(subset=['replicate_analyte_id'])
			replicate_id_df = replicate_id_df[replicate_id_df.replicate_analyte_id != None]
			replicate_id_df = replicate_id_df.sort_values(by=['replicate_analyte_id'])

			replicate_analyte_array=[]





			for analyte_id in total_df.analyte_id.unique():
				replicate_analyte_array.append(analyte_id)
		else:
			analyte_array = []









		self.ms1_rightButton = QToolButton(self)
		self.ms1_rightButton.setArrowType(Qt.RightArrow)
		self.ms1_rightButton.clicked.connect(self.is_right_button_clicked)
		self.ms1_rightButton.clicked.connect(self.Plot_ms1)


		self.ms1_leftButton = QToolButton(self)
		self.ms1_leftButton.setArrowType(Qt.LeftArrow)
		self.ms1_leftButton.clicked.connect(self.is_left_button_clicked)
		self.ms1_leftButton.clicked.connect(self.Plot_ms1)


		self.lock_aspect = QPushButton(self)

		self.lock_aspect.setCheckable(True) 
		self.lock_aspect.setIcon(QIcon('Unlock_icon.png'))
		self.lock_aspect.clicked.connect(self.Plot_ms1)
		self.lock_aspect.clicked.connect(self.Plot_ms2)




		

		# if analyteID == 'Analyte ID':


		# else:
		# 	pass	
		# if analyteID == 'replicate_analyte_id':
		# 	n = 0
		# 	while n < len(replicate_analyte_array):

		# 		self.comboAnalyte.addItem(str(replicate_analyte_array[n]))
				
		# 		n+=1	
		# else:
		# 	pass
		

		self.massRangeLabel = QLabel("m/z Range", self)

  
		self.mMinBox = QLineEdit(self)
		self.mMinBox.setText('0')
		self.mMinBox.setValidator(QDoubleValidator())

		self.mMinBox.setAlignment(Qt.AlignRight)
		# self.mMinBox.setValidator(QRegExpValidator())


		self.mMinBox.textChanged.connect(self.mMinLimit)

		if open_state == 2:
			print('===== Test Total DF')
			total_df_test = import_dataframe(output_path, input_type, comboText)


			max_mass = total_df_test[total_df_test.mz == total_df_test.mz.max()].mz
			
			max_mass_value = max_mass.to_numpy()

			max_rt = total_df_test[total_df_test.rt == total_df_test.rt.max()].rt
			
			max_rt_value = max_rt.to_numpy()



		self.mMaxBox = QLineEdit(self)


		if open_state == 2:

			self.mMaxBox.setText(str(round(max_mass_value[0])+1))

		else:
			self.mMaxBox.setText(str(0))




		self.mMaxBox.setValidator(QDoubleValidator())
		self.mMaxBox.setAlignment(Qt.AlignRight)
		self.mMaxBox.textChanged.connect(self.mMaxLimit)



		self.RTLabel = QLabel("rt Range", self)


		self.rtMinBox = QLineEdit(self)
		self.rtMinBox.setText('0')

		self.rtMinBox.setValidator(QDoubleValidator())
		self.rtMinBox.setAlignment(Qt.AlignRight)
		self.rtMinBox.textChanged.connect(self.rtMinLimit)


		self.rtMaxBox = QLineEdit(self)

		if open_state == 2:

			self.rtMaxBox.setText(str(round(max_rt_value[0])+1))

		else:
			self.rtMaxBox.setText(str(0))




		self.rtMaxBox.setValidator(QDoubleValidator())
		self.rtMaxBox.setAlignment(Qt.AlignRight)
		self.rtMaxBox.textChanged.connect(self.rtMaxLimit)


		# self.mtoLabel = QLabel("to",self)
		# self.mtoLabel.resize(50,50)
		# self.mtoLabel.move(2235,650 + move_x)

		# self.rttoLabel = QLabel("to",self)
		# self.rttoLabel.resize(50,50)
		# self.rttoLabel.move(2235,800 + move_x)

		self.massCheckbox = QCheckBox('Link m/z and rt',self)



		if self.massCheckbox.isChecked() == False:
			self.mMinBox.setText('0')
			if open_state ==2:
				self.mMaxBox.setText(str(round(max_mass_value[0])+1))


		
		# self.rtCheckbox = QCheckBox('Link m/z',self)
		# self.rtCheckbox.resize(500,50)
		# self.rtCheckbox.move(2420,800 + move_x)



		self.resetButton = QPushButton("Reset All", self)
		self.resetButton.clicked.connect(self.reset_all)








		self.resetMassButton = QPushButton("Reset", self)
		self.resetMassButton.clicked.connect(self.reset_massRange)







		self.resetRTButton = QPushButton("Reset", self)
		self.resetRTButton.clicked.connect(self.reset_rtRange)




		# if Sample_State == False:
		self.massCheckbox.stateChanged.connect(self.Plot_Sample)
		self.massCheckbox.stateChanged.connect(self.Sample_Plot_DF)
		# self.rtCheckbox.stateChanged.connect(self.Plot_Sample)
		# self.rtCheckbox.stateChanged.connect(self.Sample_Plot_DF)
		self.resetButton.clicked.connect(self.Sample_Plot_DF)
		self.resetButton.clicked.connect(self.Plot_Sample)
		self.resetMassButton.clicked.connect(self.Sample_Plot_DF)
		self.resetMassButton.clicked.connect(self.Plot_Sample)
		self.resetRTButton.clicked.connect(self.Sample_Plot_DF)
		self.resetRTButton.clicked.connect(self.Plot_Sample)
		self.comboAnalyte.currentTextChanged.connect(self.Sample_Plot_DF)
		self.comboAnalyte.currentTextChanged.connect(self.Plot_Sample)

		self.comboAnalyte.currentTextChanged.connect(self.Plot_ms1)
		self.comboAnalyte.currentTextChanged.connect(self.Plot_ms2)
		# else:
		# 	pass
		#if Sample_State == True and Replicate_State == False:
	# if Replicate_State == False:
		self.massCheckbox.stateChanged.connect(self.Plot_Replicate)
		self.massCheckbox.stateChanged.connect(self.Replicate_Plot_DF)
		# self.rtCheckbox.stateChanged.connect(self.Plot_Replicate)
		# self.rtCheckbox.stateChanged.connect(self.Replicate_Plot_DF)
		self.resetButton.clicked.connect(self.Replicate_Plot_DF)
		self.resetButton.clicked.connect(self.Plot_Replicate)
		self.resetMassButton.clicked.connect(self.Replicate_Plot_DF)
		self.resetMassButton.clicked.connect(self.Plot_Replicate)
		self.resetRTButton.clicked.connect(self.Replicate_Plot_DF)
		self.resetRTButton.clicked.connect(self.Plot_Replicate)
		self.resetRTButton.clicked.connect(self.Replicate_Plot_DF)
		self.resetRTButton.clicked.connect(self.Plot_Replicate)
		self.comboAnalyte.currentTextChanged.connect(self.Replicate_Plot_DF)
		self.comboAnalyte.currentTextChanged.connect(self.Plot_Replicate)
		# self.comboAnalyte.currentTextChanged.connect(self.Experiment_Plot_DF)
		self.comboAnalyte.currentTextChanged.connect(self.Diversity_Plot_DF)
		self.comboAnalyte.currentTextChanged.connect(self.Plot_Experiment)
		self.comboAnalyte.currentTextChanged.connect(self.Diversity_Plot_DF)
		self.comboAnalyte.currentTextChanged.connect(self.Plot_Diversity)
		# else:
		# 	pass



		self.Toggle_rt = QPushButton("Toggle Scatter Plot",self)
		self.Toggle_rt.setCheckable(True)
		self.Toggle_rt.clicked.connect(self.toggle_rt)
		self.Toggle_rt.clicked.connect(self.Plot_Sample)
		self.Toggle_rt.clicked.connect(self.Plot_Replicate)
		self.Toggle_rt.setStyleSheet("background-color : lightgrey")



		self.export_csv = QPushButton("CSV",self)
		self.export_csv.clicked.connect(self.sample_export)
		self.export_mgf = QPushButton("MGF",self)
		self.export_mgf.clicked.connect(self.ms1_DF_Export)
		# self.export_mgf = QPushButton("mgf",self)
		self.export_label_1 = QLabel('m/z and Retention Time')
		self.export_label_1.setFont(QFont('Times', 10))
		self.export_label_2 = QLabel('MS<sup>1</sup> and MS<sup>2</sup>')
		self.export_label_2.setFont(QFont('Times',10))
		self.space_label1 = QLabel('')



################ Center Left of interface (Region 3) ###################
		


		tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
		Sample_State = tab_state.return_sample_state()
		Replicate_State = tab_state.return_replicate_state()
		Experiment_State = tab_state.return_experiment_state()


		
############## Null On Off ####################

		# creating a push button 
		self.Nullbutton = QPushButton("On", self) 

		# setting geometry of button 


		# setting checkable to true 
		self.Nullbutton.setCheckable(True) 

		# setting calling method by button 

		self.Nullbutton.clicked.connect(self.nullButton)


		# setting default color of button to light-grey 
		self.Nullbutton.setStyleSheet("background-color : lightgreen") 

############## Blank On Off ####################

		# creating a push button 
		self.Blankbutton = QPushButton("On", self) 

		# setting geometry of button 


		# setting checkable to true 
		self.Blankbutton.setCheckable(True) 

		# setting calling method by button 

		self.Blankbutton.clicked.connect(self.blankButton)


		# setting default color of button to light-grey 
		self.Blankbutton.setStyleSheet("background-color : lightgreen") 

		self.Blankbutton.clicked.connect(self.Plot_Sample)
		self.Nullbutton.clicked.connect(self.Plot_Sample)

		self.Nullbutton.clicked.connect(self.Plot_Replicate)
		self.Nullbutton.clicked.connect(self.Plot_Sample)
		
		self.Nullbutton.clicked.connect(self.Plot_Replicate)
		self.Blankbutton.clicked.connect(self.Plot_Replicate)

		self.Blankbutton.clicked.connect(self.Plot_Experiment)

		self.Blankbutton.clicked.connect(self.Plot_Diversity)

		# namespace = {'pg': pg, 'np': np}
		# self.c = pyqtgraph.console.ConsoleWidget(namespace=namespace, text='test')
		# self.centralwidget = self.Window

		self.c = QtWidgets.QTextBrowser()

########### Grid Region 2 ###########################
		self.Filters = QLabel("Filter Options")
		self.Filters.setFont(QFont('Times', 12))
		self.Export = QLabel("Export Options")
		self.Export.setFont(QFont('Times', 12))
		self.Line1 = QLabel("__________________________________________________")
		self.Line2 = QLabel("__________________________________________________")
		self.Line3 = QLabel("__________________________________________________")
		self.Line4 = QLabel("__________________________________________________")
		self.Line5 = QLabel("__________________________________________________")
		self.Line6 = QLabel("__________________________________________________")
		self.Blank1 = QLabel("   ")
		layoutRegion2 = QGridLayout()

		# layoutRegion2.setRowStretch(0, 1)
		layoutRegion2.addWidget(self.chooseSample,0,0)
		layoutRegion2.addWidget(self.comboBox,0,1)
		layoutRegion2.addWidget(self.chooseAnalyte,1,0)
		layoutRegion2.addWidget(self.comboAnalyte,1,1)

		layoutRegion4 = QGridLayout()
		layoutRegion4.addWidget(self.Blank1,0,0,1,3,QtCore.Qt.AlignCenter)
		layoutRegion4.addWidget(self.Filters,1,0,1,3)
		layoutRegion4.addWidget(self.Line1,2,0,1,3,QtCore.Qt.AlignTop)
		layoutRegion4.addWidget(self.showNull,3,0)
		layoutRegion4.addWidget(self.Nullbutton,3,1)
		layoutRegion4.addWidget(self.showBlanks,4,0)
		layoutRegion4.addWidget(self.Blankbutton,4,1)


		layoutRegion4.addWidget(self.Line2,5,0,1,3)
		layoutRegion4.addWidget(self.massRangeLabel,6,0)
		layoutRegion4.addWidget(self.resetMassButton,6,1)
		layoutRegion4.addWidget(self.mMinBox,7,0)
		layoutRegion4.addWidget(self.mMaxBox,7,1)


		layoutRegion4.addWidget(self.Line3,8,0,1,3)
		layoutRegion4.addWidget(self.RTLabel,9,0)
		layoutRegion4.addWidget(self.resetRTButton,9,1)
		layoutRegion4.addWidget(self.rtMinBox,10,0)
		layoutRegion4.addWidget(self.rtMaxBox,10,1)
		layoutRegion4.addWidget(self.massCheckbox,12,1)

		layoutRegion4.addWidget(self.Line4,11,0,1,3)
		layoutRegion4.addWidget(self.resetButton,13,0)
		layoutRegion4.addWidget(self.Toggle_rt,12,0)
		layoutRegion4.addWidget(self.Blank1,14,0,1,3)
		layoutRegion4.addWidget(self.Export,15,0,1,3)
		layoutRegion4.addWidget(self.Line5,16,0,1,3)
		layoutRegion4.addWidget(self.export_label_1,17,0,QtCore.Qt.AlignCenter)
		layoutRegion4.addWidget(self.export_csv,18,0)
		layoutRegion4.addWidget(self.export_label_2,17,1,QtCore.Qt.AlignCenter)
		layoutRegion4.addWidget(self.export_mgf,18,1)
		layoutRegion4.addWidget(self.Line6,19,0,1,3)
		layoutRegion4.addWidget(self.c,20,0,1,2)
		# layoutRegion4.addWidget(self.space_label1,0,0)





		comboText = self.comboBox.currentText()


################## Sample ###########################

		self.mz_label1 = pg.LabelItem(justify='right')
		self.mz_plot = pg.PlotWidget()
		
		self.mz_plot.setLabel(axis='bottom',text='m/z',**{ 'font-size': '10pt'}, **{'font': 'Italicized'})
		self.mz_plot.setLabel(axis='bottom',text='m/z',**{ 'font-size': '10pt'}, **{'font': 'Italicized'})
		
		self.mz_plot.setLabel(axis='left',text='Intensity',**{ 'font-size': '10pt'})
		self.mz_plot.setMouseEnabled(x=True, y=False)
		self.mz_plot.setStyleSheet("background-color: white;  border:1px solid black")
		self.mz_plot.setTitle("No Data Loaded")




		self.rt_plot = pg.PlotWidget()
		self.rt_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		self.rt_plot.setLabel(axis='bottom',text='rt',**{ 'font-size': '10pt'})
		self.rt_plot.setLabel(axis='left',text='Intensity',**{ 'font-size': '10pt'})
		
		self.rt_plot.setMouseEnabled(x=True, y=False)


		self.rt_plot.setTitle("No Data Loaded")
		self.ms1_plot = pg.PlotWidget()
		self.ms1_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		self.ms1_plot.setLabel(axis='bottom',text='m/z',**{ 'font-size': '10pt'})
		self.ms1_plot.setTitle("MS<sup>1</sup> Spectra")
		self.ms1_plot.setLabel(axis='left',text='Intensity',**{ 'font-size': '10pt'})
		self.ms1_plot.setMouseEnabled(x=True, y=False)
		self.ms2_plot = pg.PlotWidget()
		self.ms2_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		self.ms2_plot.setLabel(axis='bottom',text='m/z',**{ 'font-size': '10pt'}, **{'font': 'Italicized'})
		self.ms2_plot.setTitle("MS<sup>2</sup> Spectra")
		self.ms2_plot.setLabel(axis='left',text='Intensity',**{ 'font-size': '10pt'})
		self.ms2_plot.setMouseEnabled(x=True, y=False)


################################################# Replicate

############# Replicate #############################


		self.mz_r1_plot = pg.PlotWidget()
		self.mz_r1_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		self.mz_r1_plot.setTitle("Replicate 1 m/z vs intensity")
		self.mz_r1_plot.setLabel(axis='bottom',text='m/z',**{ 'font-size': '10pt'}, **{'font': 'Italicized'})
		self.mz_r1_plot.setLabel(axis='left',text='Intensity')
		self.mz_r1_plot.setMouseEnabled(x=True, y=False)
		self.mz_r2_plot = pg.PlotWidget()
		self.mz_r2_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		self.mz_r2_plot.setTitle("Replicate 2 m/z vs intensity")
		self.mz_r2_plot.setLabel(axis='bottom',text='m/z',**{ 'font-size': '10pt'}, **{'font': 'Italicized'})
		self.mz_r2_plot.setLabel(axis='left',text='Intensity')
		self.mz_r2_plot.setMouseEnabled(x=True, y=False)
		self.mz_r3_plot = pg.PlotWidget()
		self.mz_r3_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		self.mz_r3_plot.setTitle("Replicate 3 m/z vs intensity")
		self.mz_r3_plot.setLabel(axis='bottom',text='m/z',**{ 'font-size': '10pt'}, **{'font': 'Italicized'})
		self.mz_r3_plot.setLabel(axis='left',text='Intensity')
		self.mz_r3_plot.setMouseEnabled(x=True, y=False)
		self.rt_r1_plot = pg.PlotWidget()
		self.rt_r1_plot.setTitle("Replicate 1 rt vs intesnity")
		self.rt_r1_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		self.rt_r1_plot.setLabel(axis='bottom',text='rt',**{ 'font-size': '10pt'})
		self.rt_r1_plot.setLabel(axis='left',text='Intensity')
		self.rt_r1_plot.setMouseEnabled(x=True, y=False)
		self.rt_r2_plot = pg.PlotWidget()
		self.rt_r2_plot.setTitle("Replicate 2 rt vs intesnity")
		self.rt_r2_plot.setLabel(axis='bottom',text='rt',**{ 'font-size': '10pt'})
		self.rt_r2_plot.setLabel(axis='left',text='Intensity')
		self.rt_r2_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		self.rt_r2_plot.setMouseEnabled(x=True, y=False)
		self.rt_r3_plot = pg.PlotWidget()
		self.rt_r3_plot.setTitle("Replicate 3 rt vs intesnity")
		self.rt_r3_plot.setLabel(axis='bottom',text='rt',**{ 'font-size': '10pt'})
		self.rt_r3_plot.setLabel(axis='left',text='Intensity')
		self.rt_r3_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		self.rt_r3_plot.setMouseEnabled(x=True, y=False)



		# self.No_Data_Label1 = QLabel('No Data Loaded')
		# self.No_Data_Label2 = QLabel('No Data Loaded')
		# self.No_Data_Label3 = QLabel('No Data Loaded')
		# self.No_Data_Label4 = QLabel('No Data Loaded')
		# self.No_Data_Label5 = QLabel('No Data Loaded')
		# self.No_Data_Label6 = QLabel('No Data Loaded')
		# self.No_Data_Label7 = QLabel('No Data Loaded')
		# self.No_Data_Label8 = QLabel('No Data Loaded')



############# Experiment and Diversity #############################

		self.ex_bOn_plot = pg.PlotWidget()
		self.ex_bOn_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")





		self.div_plot = pg.PlotWidget()
		self.div_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		# self.ex_bOn_plot.move(20,250 + move_y)

		if open_state ==2:
			if len(sample_names) < 16:
				self.ex_bOn_plot.resize(2000,100*len(sample_names))
			else:
				self.ex_bOn_plot.resize(2000,1500)
				self.ex_bOn_plot.setRange(yRange=(0,15), padding = 0)
			self.ex_bOn_plot.setLabel(axis='bottom',text='Max Mass')



			if len(sample_names) < 16:
				self.div_plot.resize(2000,100*len(sample_names))
			else:
				self.div_plot.resize(2000,1500)
				self.div_plot.setRange(yRange=(0,15), padding = 0)
			self.div_plot.setLabel(axis='bottom',text='Analyte ID')
##############################################################


		self.analyte_legend = pg.PlotWidget(self)
		self.analyte_legend.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		self.analyte_legend.hideAxis('left')
		self.analyte_legend.hideAxis('bottom')
		self.analyte_legend.setMouseEnabled(x=True, y=True)
		layoutRegion9 = QGridLayout()

		layoutRegion9.addWidget(self.analyte_legend,0,1)
		# layoutRegion9.addWidget(self.s1,0,1)

##############################################################
		self.ms2analytelabel = QLabel(self)
		pixmap = QPixmap('MS2Analyte_logo.png')

		self.ms2analytelabel.setPixmap(pixmap)
		layoutRegion7 = QGridLayout()

		layoutRegion7.addWidget(self.ms2analytelabel)
##############################################################

############ Grid Region 3 ###################################

		layoutRegion3 = QGridLayout()
		layoutRegion3.addWidget(self.mz_plot,0,1)
		layoutRegion3.addWidget(self.rt_plot,1,1)
		layoutRegion3.addWidget(self.mz_r1_plot,1,1)
		layoutRegion3.addWidget(self.mz_r2_plot,2,1)
		layoutRegion3.addWidget(self.mz_r3_plot,3,1)
		layoutRegion3.addWidget(self.rt_r1_plot,4,1)
		layoutRegion3.addWidget(self.rt_r2_plot,5,1)
		layoutRegion3.addWidget(self.rt_r3_plot,6,1)

		# layoutRegion3.addWidget(self.No_Data_Label1,0,1)
		# layoutRegion3.addWidget(self.No_Data_Label2,1,1)
		# layoutRegion3.addWidget(self.No_Data_Label3,1,1)
		# layoutRegion3.addWidget(self.No_Data_Label4,2,1)
		# layoutRegion3.addWidget(self.No_Data_Label5,3,1)
		# layoutRegion3.addWidget(self.No_Data_Label6,4,1)
		# layoutRegion3.addWidget(self.No_Data_Label7,5,1)
		# layoutRegion3.addWidget(self.No_Data_Label8,6,1)

		layoutRegion5 = QGridLayout()
		layoutRegion5.addWidget(self.ms1_plot,0,1)
		layoutRegion5.addWidget(self.ms2_plot,1,1)








		layoutRegion6 = QGridLayout()

		layoutRegion6.addWidget(self.ex_bOn_plot,1,0)

		layoutRegion8 = QGridLayout()

		layoutRegion8.addWidget(self.div_plot,1,0)



		# layoutRegion10 = QGridLayout()
		# layoutRegion10.addWidget(self.ms1_leftButton,2,1,QtCore.Qt.AlignLeft)
		# layoutRegion10.addWidget(self.ms1_rightButton,2,1,QtCore.Qt.AlignRight)

		layoutRegion11 = QGridLayout()
		# layoutRegion11.addWidget(self.lock_aspect,2,2,QtCore.Qt.AlignCenter)
		# layoutRegion11.addWidget(self.ms1_rightButton,2,1,QtCore.Qt.AlignCenter)
		# layoutRegion11.addWidget(self.ms1_leftButton,2,0,QtCore.Qt.AlignCenter)

		layoutRegion11.addWidget(self.lock_aspect,2,2,QtCore.Qt.AlignTop)
		layoutRegion11.addWidget(self.ms1_rightButton,2,1,QtCore.Qt.AlignTop)
		layoutRegion11.addWidget(self.ms1_leftButton,2,0,QtCore.Qt.AlignTop)
		layoutRegion11.setContentsMargins(0,5, 5, 0)

		# layoutRegion12 = QGridLayout()
		# layoutRegion12.addWidget(self.menu_bar,0,0)


############# Create Main Grid ###############################

		# g.addWidget(button, N, 0, 1, N, QtCore.Qt.AlignCenter)
		layoutMain.setColumnStretch(0, 4)
		layoutMain.setColumnStretch(1, 2)
		# layoutMain.setRowStretch(1, 0.5)
		# layoutRegion2.setColumnStretch(0,4)
		layoutMain.addLayout(layoutRegion3,2,0,QtCore.Qt.AlignRight)
		layoutMain.addLayout(layoutRegion5,2,1,QtCore.Qt.AlignRight)
		layoutMain.addLayout(layoutRegion1,1,0,QtCore.Qt.AlignLeft)

		# layoutMain.setSpacing(0)
		layoutMain.addLayout(layoutRegion2,1,1, 1, 1, QtCore.Qt.AlignRight)
		layoutMain.addLayout(layoutRegion7,1,3, QtCore.Qt.AlignCenter)
		layoutMain.addLayout(layoutRegion4,2,3, 1, 1,QtCore.Qt.AlignTop)
		layoutMain.addLayout(layoutRegion6,2,0, 1,1)
		layoutMain.addLayout(layoutRegion8,2,0, 1,1)

		layoutMain.addLayout(layoutRegion9,2,2)
		# layoutMain.addLayout(layoutRegion10,1,1,QtCore.Qt.AlignTop)
		layoutMain.addLayout(layoutRegion11,2,1,1, 1, QtCore.Qt.AlignRight)
		# layoutMain.addLayout(layoutRegion12,0,0)

############################################################################



		# layoutMain.addLayout(layoutRegion3_R,1,0)
		# layoutMain.addLayout(layoutRegion3_E,1,0)
		self.horizontalGroupBox.setLayout(layoutMain)





		self.ms1_rightButton.setEnabled(False)
		self.ms1_leftButton.setEnabled(False)

################### Plots #################################


	def loading_sample(self):
		self.c.clear()
		self.c.append('Loading Sample Tab...')


	def loading_replicate(self):
		self.c.clear()
		self.c.append('Loading Replicate Tab...')

	def loading_experiment(self):
		self.c.clear()
		self.c.append('Loading Experiment Tab...')

	def loading_diversity(self):
		self.c.clear()
		self.c.append('Loading Diversity Tab...')

	def fill_combobox(self):
		if open_state == 2:
			comboText = self.comboBox.currentText()
			chooseAnalyte = self.chooseAnalyte.currentText()
			print('FILL COMBO')
			print(comboText)
			total_df = import_dataframe(output_path, input_type, comboText)
			current_sample = Current_Sample(comboText)

			
			# sorted_df = total_df['analyte_id'].replace('',np.nan, inplace=True)
			# sorted_df = total_df.dropna(subset=['analyte_id'],inplace=True)
			# print(sorted_df)


			replicate_id_df = total_df.dropna(subset=['replicate_analyte_id'])
			experiment_id_df = total_df.dropna(subset=['experiment_analyte_id'])
			sorted_df = total_df.dropna(subset=['analyte_id'])
			sorted_df=total_df[total_df.analyte_id != None]
			sorted_df = sorted_df.sort_values(by=['analyte_id'])
			# total_df['analyte_id'].replace('',np.nan, inplace=True)
			mz_array=[]
			analyte_array=[]

			for analyte_id in sorted_df.analyte_id.unique():
				mz_array.append([sorted_df[sorted_df["analyte_id"] == analyte_id]["mz"].to_numpy(),
			                     sorted_df[sorted_df["analyte_id"] == analyte_id]["intensity"].to_numpy()])

				analyte_array.append(analyte_id)


			replicate_id_df = replicate_id_df.dropna(subset=['replicate_analyte_id'])
			replicate_id_df = replicate_id_df[replicate_id_df.replicate_analyte_id != None]
			replicate_id_df = replicate_id_df.sort_values(by=['replicate_analyte_id'])

			replicate_analyte_array=[]



			for analyte_id in replicate_id_df.replicate_analyte_id.unique():
				replicate_analyte_array.append(analyte_id)


			experiment_id_df = experiment_id_df.dropna(subset=['experiment_analyte_id'])
			experiment_id_df = experiment_id_df[experiment_id_df.experiment_analyte_id != None]
			experiment_id_df = experiment_id_df.sort_values(by=['experiment_analyte_id'])

			experiment_analyte_array=[]



			for analyte_id in experiment_id_df.experiment_analyte_id.unique():
				experiment_analyte_array.append(analyte_id)
			# self.comboAnalyte.clear()



			if chooseAnalyte == "Sample Analyte ID":


				self.comboAnalyte.blockSignals(True)
				self.comboAnalyte.clear()
				self.comboAnalyte.addItem("Show All")
				n = 0
				while n < len(analyte_array):

					self.comboAnalyte.addItem(str(analyte_array[n]))
					
					n+=1
					self.comboAnalyte.blockSignals(False)
			else:
				pass



			if chooseAnalyte == "Replicate Analyte ID":


				self.comboAnalyte.blockSignals(True)
				self.comboAnalyte.clear()
				self.comboAnalyte.addItem("Show All")

				n = 0
				while n < len(replicate_analyte_array):

					self.comboAnalyte.addItem(str(replicate_analyte_array[n]))
					
					n+=1
					self.comboAnalyte.blockSignals(False)
			
			else:
				pass
			if chooseAnalyte == "Experiment Analyte ID":


				self.comboAnalyte.blockSignals(True)
				self.comboAnalyte.clear()
				self.comboAnalyte.addItem("Show All")

				n = 0
				while n < len(experiment_analyte_array):

					self.comboAnalyte.addItem(str(experiment_analyte_array[n]))
					
					n+=1
					self.comboAnalyte.blockSignals(False)
			
			else:
				pass
	

	############################## Clicking Data Points On Plots
	
	def clicked(self, plot, points):


		global lastClicked




		for p in lastClicked:


			p.resetPen()



		self.c.append('_____________________________________________')
		self.c.append('\nData Selected')

		comboText=self.comboBox.currentText()
		print('TestTest Clicked')
		total_df = import_dataframe(output_path, input_type, comboText)


		indexes = []

		for p in points:

			p.resetPen()

			p.setPen(clickedPen)

			p = p.pos()

			x = float(p[0])

			y = float(p[1])

			info_df = total_df[total_df.mz == x]
			info_df = info_df[info_df.intensity == y]

			analyte_id = info_df.analyte_id.to_numpy()

			retention_time = info_df.rt.to_numpy()

			drift_time = info_df.drift.to_numpy()

			# blank = info_df.blank_analyte.to_numpy()
			# print('Clicked Test 20')
			
			
	
			self.c.append('---------------------------------')
			self.c.append('Analyte ID: '+ str(analyte_id[0]))
			self.c.append('---------------------------------')
			self.c.append('m/z: '+str(x)+' \nIntensity: '+str(y)+' \nRetention Time: '+str(retention_time[0])+' \nDrift Time: '+str(drift_time[0]))
			self.mz_plot.setToolTip('Data Selected\n'+'\nAnalyte ID: '+ str(analyte_id[0])+'\nm/z: '+str(x)+' \nIntensity: '+str(y)+' \nRetention Time: '+str(retention_time[0])+' \nDrift Time: '+str(drift_time[0]))

			self.mz_r1_plot.setToolTip('Data Selected\n'+'\nAnalyte ID: '+ str(analyte_id[0])+'\nm/z: '+str(x)+' \nIntensity: '+str(y)+' \nRetention Time: '+str(retention_time[0])+' \nDrift Time: '+str(drift_time[0]))
			self.mz_r2_plot.setToolTip('Data Selected\n'+'\nAnalyte ID: '+ str(analyte_id[0])+'\nm/z: '+str(x)+' \nIntensity: '+str(y)+' \nRetention Time: '+str(retention_time[0])+' \nDrift Time: '+str(drift_time[0]))
			self.mz_r3_plot.setToolTip('Data Selected\n'+'\nAnalyte ID: '+ str(analyte_id[0])+'\nm/z: '+str(x)+' \nIntensity: '+str(y)+' \nRetention Time: '+str(retention_time[0])+' \nDrift Time: '+str(drift_time[0]))

			# text = pg.TextItem('Data Selected\n'+'\nAnalyte ID: '+ str(analyte_id[0])+'\nm/z: '+str(x)+' \nIntensity: '+str(y)+' \nRetention Time: '+str(retention_time[0])+' \nDrift Time: '+str(drift_time[0])+' \nBlank: '+str(blank[0]),color=(0,0,0))
			# self.mz_plot.addItem(text)


			# text.setPos(1000, 300000)
		lastClicked = points

	def clicked_Legend(self, plot, points):


		global lastClicked_legend




		for p in lastClicked_legend:


			p.resetPen()





		comboText=self.comboBox.currentText()




		indexes = []

		for p in points:

			p.resetPen()

			p.setPen(clickedPen_legend)

			p = p.pos()

			x = float(p[0])

			y = float(p[1])


			# blank = info_df.blank_analyte.to_numpy()
			# print('Clicked Test 20')
			
			
	


			# text = pg.TextItem('Data Selected\n'+'\nAnalyte ID: '+ str(analyte_id[0])+'\nm/z: '+str(x)+' \nIntensity: '+str(y)+' \nRetention Time: '+str(retention_time[0])+' \nDrift Time: '+str(drift_time[0])+' \nBlank: '+str(blank[0]),color=(0,0,0))
			# self.mz_plot.addItem(text)


			# text.setPos(1000, 300000)
		lastClicked_legend = points


	def ms1_clicked(self, plot, points):
		indexes = []

		for p in points:
			print("ms1 CLICKED")

	def Plot_Sample(self):
		if open_state == 2:


			print('OPEN STATE',open_state)
			global lastClicked
			lastClicked = []
			print('PLOT SAMPLE BEGIN')
			
			print('TEST11111111')
			tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
			print('TEST2222222')
			Sample_State = tab_state.return_sample_state()
			print('TEST3333333')
			Replicate_State = tab_state.return_replicate_state()
			print('TEST44444444')
			Experiment_State = tab_state.return_experiment_state()
			print("TESTTESTTEST")
			print('Loading Sample Tab**')
			state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), 0, 1200, 0, 7,self.Toggle_rt.isChecked())
			Toggle_rt = state.return_Toggle_rt()

			if Sample_State == False:

				self.c.update()
				chooseAnalyte = self.chooseAnalyte.currentText()
				self.mz_plot.disableAutoRange(axis=None)
				self.rt_plot.disableAutoRange(axis=None)
				comboAnalyte=self.comboAnalyte.currentText()
				comboText = self.comboBox.currentText()
				self.mz_plot.setTitle(comboText)
				self.rt_plot.setTitle(comboText)

				print('Loading Sample Tab***')
				
				############### Fetch data for mz_plot and plot sample from comboBox ###############

				sorted_df_mz, sorted_df_rt= self.Sample_Plot_DF()
				print('ANALYTE ID 2',sorted_df_mz[sorted_df_mz.analyte_id == 2])



				mz_array = []
				analyte_mz_array = []
				rt_array = []
				analyte_rt_array = []









					# df=df[df.scan == max_scan_value[0]]
					# df.sort_values('scan', inplace=True)

				if chooseAnalyte == "Sample Analyte ID":
					for analyte_id in sorted_df_mz.analyte_id.unique():

						mz_array.append([sorted_df_mz[sorted_df_mz["analyte_id"] == analyte_id]["mz"].to_numpy(),
					                     sorted_df_mz[sorted_df_mz["analyte_id"] == analyte_id]["intensity"].to_numpy()])

						analyte_mz_array.append(analyte_id)



					if Toggle_rt == True:
					

						for analyte_id in sorted_df_rt.analyte_id.unique():
							rt_array.append([sorted_df_rt[sorted_df_rt["analyte_id"] == analyte_id]["rt"].to_numpy(),
						                     sorted_df_rt[sorted_df_rt["analyte_id"] == analyte_id]["intensity"].to_numpy()])
							analyte_rt_array.append(analyte_id)
					else:



						for analyte_id in sorted_df_rt.analyte_id.unique():

							analyte_df = sorted_df_rt[sorted_df_rt.analyte_id == analyte_id]

							max_scan_df = analyte_df[analyte_df.intensity == analyte_df.intensity.max()]

							max_peak_id_array = max_scan_df.peak_id.to_numpy()
							max_peak_id = int(max_peak_id_array[0])




							sorted_df_rt1 = analyte_df[analyte_df.peak_id == int(max_peak_id)]


							rt_array.append([sorted_df_rt1[sorted_df_rt1["analyte_id"] == analyte_id]["rt"].to_numpy(),
						                     sorted_df_rt1[sorted_df_rt1["analyte_id"] == analyte_id]["intensity"].to_numpy()])
							analyte_rt_array.append(analyte_id)


					analyte_mz_array = [0 if str(i) == 'nan' else i for i in analyte_mz_array]

					analyte_rt_array = [0 if str(i) == 'nan' else i for i in analyte_rt_array]



					analyte_colour_mz_array = [0]*len(analyte_mz_array)
					analyte_colour_rt_array = [0]*len(analyte_rt_array)

					i = 0
					while i < len(analyte_mz_array):
						analyte_colour_mz_array[i] = int((analyte_mz_array[i] + 1) % 299)
						i+=1


					i = 0
					while i < len(analyte_rt_array):

						analyte_colour_rt_array[i] = int((analyte_rt_array[i] + 1) % 299)
						i+=1
					print('Loading Sample Tab******')
				if chooseAnalyte == "Replicate Analyte ID":

					for analyte_id in sorted_df_mz.replicate_analyte_id.unique():

						mz_array.append([sorted_df_mz[sorted_df_mz["replicate_analyte_id"] == analyte_id]["mz"].to_numpy(),
					                     sorted_df_mz[sorted_df_mz["replicate_analyte_id"] == analyte_id]["intensity"].to_numpy()])

						analyte_mz_array.append(analyte_id)



					if Toggle_rt == True:
					

						for analyte_id in sorted_df_rt.analyte_id.unique():
							rt_array.append([sorted_df_rt[sorted_df_rt["analyte_id"] == analyte_id]["rt"].to_numpy(),
						                     sorted_df_rt[sorted_df_rt["analyte_id"] == analyte_id]["intensity"].to_numpy()])
							analyte_rt_array.append(analyte_id)
					else:



						for analyte_id in sorted_df_rt.analyte_id.unique():

							analyte_df = sorted_df_rt[sorted_df_rt.analyte_id == analyte_id]

							max_scan_df = analyte_df[analyte_df.intensity == analyte_df.intensity.max()]

							max_peak_id_array = max_scan_df.peak_id.to_numpy()
							max_peak_id = int(max_peak_id_array[0])



							sorted_df_rt1 = analyte_df[analyte_df.peak_id == int(max_peak_id)]


							rt_array.append([sorted_df_rt1[sorted_df_rt1["analyte_id"] == analyte_id]["rt"].to_numpy(),
						                     sorted_df_rt1[sorted_df_rt1["analyte_id"] == analyte_id]["intensity"].to_numpy()])
							analyte_rt_array.append(analyte_id)


					analyte_mz_array = [0 if str(i) == 'nan' else i for i in analyte_mz_array]

					analyte_rt_array = [0 if str(i) == 'nan' else i for i in analyte_rt_array]

					analyte_mz_array = [0 if str(i) == 'None' else i for i in analyte_mz_array]

					analyte_rt_array = [0 if str(i) == 'None' else i for i in analyte_rt_array]

					analyte_colour_mz_array = [0]*len(analyte_mz_array)
					analyte_colour_rt_array = [0]*len(analyte_rt_array)


					i = 0
					while i < len(analyte_mz_array):
						analyte_colour_mz_array[i] = int((analyte_mz_array[i] + 1) % 299)

						i+=1


					i = 0
					while i < len(analyte_rt_array):

						analyte_colour_rt_array[i] = int((analyte_rt_array[i] + 1) % 299)
						i+=1

					print('Loading Sample Tab')

				if chooseAnalyte == "Experiment Analyte ID":

					for analyte_id in sorted_df_mz.experiment_analyte_id.unique():

						mz_array.append([sorted_df_mz[sorted_df_mz["experiment_analyte_id"] == analyte_id]["mz"].to_numpy(),
					                     sorted_df_mz[sorted_df_mz["experiment_analyte_id"] == analyte_id]["intensity"].to_numpy()])

						analyte_mz_array.append(analyte_id)



					if Toggle_rt == True:
					

						for analyte_id in sorted_df_rt.analyte_id.unique():
							rt_array.append([sorted_df_rt[sorted_df_rt["analyte_id"] == analyte_id]["rt"].to_numpy(),
						                     sorted_df_rt[sorted_df_rt["analyte_id"] == analyte_id]["intensity"].to_numpy()])
							analyte_rt_array.append(analyte_id)
					else:



						for analyte_id in sorted_df_rt.analyte_id.unique():

							analyte_df = sorted_df_rt[sorted_df_rt.analyte_id == analyte_id]

							max_scan_df = analyte_df[analyte_df.intensity == analyte_df.intensity.max()]

							max_peak_id_array = max_scan_df.peak_id.to_numpy()
							max_peak_id = int(max_peak_id_array[0])


							# print(analyte_df.peak_id)


							sorted_df_rt1 = analyte_df[analyte_df.peak_id == int(max_peak_id)]


							rt_array.append([sorted_df_rt1[sorted_df_rt1["analyte_id"] == analyte_id]["rt"].to_numpy(),
						                     sorted_df_rt1[sorted_df_rt1["analyte_id"] == analyte_id]["intensity"].to_numpy()])
							analyte_rt_array.append(analyte_id)


					analyte_mz_array = [0 if str(i) == 'nan' else i for i in analyte_mz_array]

					analyte_rt_array = [0 if str(i) == 'nan' else i for i in analyte_rt_array]

					analyte_mz_array = [0 if str(i) == 'None' else i for i in analyte_mz_array]

					analyte_rt_array = [0 if str(i) == 'None' else i for i in analyte_rt_array]

					analyte_colour_mz_array = [0]*len(analyte_mz_array)
					analyte_colour_rt_array = [0]*len(analyte_rt_array)


					i = 0
					while i < len(analyte_mz_array):
						analyte_colour_mz_array[i] = int((analyte_mz_array[i] + 1) % 299)

						i+=1


					i = 0
					while i < len(analyte_rt_array):

						analyte_colour_rt_array[i] = int((analyte_rt_array[i] + 1) % 299)
						i+=1

					print('Loading Sample Tab')




				self.rt_plot.clear()
				self.mz_plot.clear()



				if Toggle_rt == False:

					if comboAnalyte != "Show All":
						alpha=1
						print(comboAnalyte)
						print('Loading Sample Tab..*')

						colour=0
						for analyte in mz_array:


							#self.mz_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_array[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_array[colour],values=5,alpha=alpha),symbol='o', symbolSize=3)
							#colour+=1

							if analyte_mz_array[colour] == float(comboAnalyte):


								self.test = self.mz_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array[colour],  values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_mz_array[colour],values=5,alpha=255),symbol='o', symbolSize=3)
								self.test.sigPointsClicked.connect(self.clicked)
							else:
								pass


							self.test = self.mz_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_mz_array[colour],values=5,alpha=alpha),symbol='o', symbolSize=3)
							self.test.sigPointsClicked.connect(self.clicked)
							colour+=1

						colour=0
						for analyte in rt_array:


							if analyte_rt_array[colour] == float(comboAnalyte):


								self.rt_plot.plot( x=analyte[0], y=analyte[1], fillLevel=-0.3,pen=pg.intColor(analyte_colour_rt_array[colour], values=5,alpha=255),brush=pg.intColor(analyte_colour_rt_array[colour], alpha=50))
								

							else:
								pass

							self.rt_plot.plot( x=analyte[0], y=analyte[1], fillLevel=-0.3,pen=pg.intColor(analyte_colour_rt_array[colour], values=5,alpha=alpha),brush=pg.intColor(analyte_colour_rt_array[colour], alpha=alpha))
							colour+=1


						print('Loading Sample Tab...')
					else:
						alpha=255
						colour=0
						for analyte in mz_array:


							self.test = self.mz_plot.plot( title="test", x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array[colour],  values=5), symbolBrush=pg.intColor( analyte_colour_mz_array[colour],values=5),symbol='o', symbolSize=3)
							#test.setAlpha(0, False)
							self.test.sigPointsClicked.connect(self.clicked)
							colour+=1
						colour=0
						for analyte in rt_array:





							
							self.rt_plot.plot( x=analyte[0], y=analyte[1], fillLevel=-0.3, pen=pg.intColor(analyte_colour_rt_array[colour], values=5),brush=pg.intColor(analyte_colour_rt_array[colour], alpha=50))
							colour+=1
						print('Loading Sample Tab...')




				else:
					if comboAnalyte != "Show All":
						alpha=1
						print('Loading Sample Tab..')

						colour=0
						for analyte in mz_array:


							#self.mz_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_array[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_array[colour],values=5,alpha=alpha),symbol='o', symbolSize=3)
							#colour+=1

							if analyte_mz_array[colour] == float(comboAnalyte):


								self.test = self.mz_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array[colour],  values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_mz_array[colour],values=5,alpha=255),symbol='o', symbolSize=3)
								self.test.sigPointsClicked.connect(self.clicked)
							else:
								pass


							self.test = self.mz_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_mz_array[colour],values=5,alpha=alpha),symbol='o', symbolSize=3)
							self.test.sigPointsClicked.connect(self.clicked)
							colour+=1

						colour=0
						for analyte in rt_array:


							if analyte_rt_array[colour] == float(comboAnalyte):


								self.rt_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array[colour],  values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_rt_array[colour], values=5, alpha=255),symbol='o', symbolSize=3)
								

							else:
								pass

							self.rt_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_rt_array[colour], values=5, alpha=alpha),symbol='o', symbolSize=3)
							colour+=1


						print('Loading Sample Tab...')
					else:
						alpha=255
						colour=0
						for analyte in mz_array:


							self.test = self.mz_plot.plot( title="test", x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array[colour],  values=5), symbolBrush=pg.intColor( analyte_colour_mz_array[colour],values=5),symbol='o', symbolSize=3)
							#test.setAlpha(0, False)
							self.test.sigPointsClicked.connect(self.clicked)
							colour+=1
						colour=0
						for analyte in rt_array:





							self.rt_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array[colour],  values=5), symbolBrush=pg.intColor( analyte_colour_rt_array[colour], values=5),symbol='o', symbolSize=3)

							colour+=1
						print('Loading Sample Tab...')
				



				#analyte_array[0] = 0
				#print(analyte_array)






				self.mz_plot.enableAutoRange(axis=None)
				self.rt_plot.enableAutoRange(axis=None)
				state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), self.rtMinBox.text(), self.rtMaxBox.text(), self.Toggle_rt.isChecked())

				self.mz_plot.setXRange(float(state.return_massMin_state()),float(state.return_massMax_state()))
				self.rt_plot.setXRange(float(state.return_rtMin_state()),float(state.return_rtMax_state()))
			else:
				pass






			
			self.c.update()
			self.c.append('Sample Tab Loaded Succesfully')
			print('Sample Tab Loaded Succesfully')
			print('PLOT SAMPLE END')
		else:
			pass


##################################################################################################################	

		return 


	def Plot_Replicate(self):





		global lastClicked
		lastClicked = []
		print('PLOT REPLICATE BEGIN')
		print('Loading Replicate Tab')
		chooseAnalyte = self.chooseAnalyte.currentText()
		tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
		Sample_State = tab_state.return_sample_state()
		Replicate_State = tab_state.return_replicate_state()
		Experiment_State = tab_state.return_experiment_state()
		state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), 0, 1200, 0, 7,self.Toggle_rt.isChecked())
		Toggle_rt = state.return_Toggle_rt()


		print('Loading Replicate Tab')
		if open_state ==2:
			if Replicate_State == False:


				self.mz_r1_plot.disableAutoRange(axis=None)
				self.mz_r2_plot.disableAutoRange(axis=None)
				self.mz_r3_plot.disableAutoRange(axis=None)
				self.rt_r1_plot.disableAutoRange(axis=None)
				self.rt_r2_plot.disableAutoRange(axis=None)
				self.rt_r3_plot.disableAutoRange(axis=None)
				comboAnalyte=self.comboAnalyte.currentText()
			

				
				############### Fetch data for mz_plot and plot sample from comboBox ###############

				sorted_df_mz_r1, sorted_df_rt_r1,sorted_df_mz_r2, sorted_df_rt_r2,sorted_df_mz_r3, sorted_df_rt_r3 = self.Replicate_Plot_DF()


				mz_array_r1 = []
				analyte_mz_array_r1 = []
				rt_array_r1 = []
				analyte_rt_array_r1 = []


				if chooseAnalyte == "Sample Analyte ID":



					print('Loading Replicate Tab for Sample Analyte ID')
					for analyte_id in sorted_df_mz_r1.analyte_id.unique():

						mz_array_r1.append([sorted_df_mz_r1[sorted_df_mz_r1["analyte_id"] == analyte_id]["mz"].to_numpy(),
					                     sorted_df_mz_r1[sorted_df_mz_r1["analyte_id"] == analyte_id]["intensity"].to_numpy(),])

						analyte_mz_array_r1.append(analyte_id)






					if Toggle_rt == True:
					
						print('HERE')
						for analyte_id in sorted_df_rt_r1.analyte_id.unique():
							rt_array_r1.append([sorted_df_rt_r1[sorted_df_rt_r1["analyte_id"] == analyte_id]["rt"].to_numpy(),
						                     sorted_df_rt_r1[sorted_df_rt_r1["analyte_id"] == analyte_id]["intensity"].to_numpy()])
							analyte_rt_array_r1.append(analyte_id)

					else:

						

						for analyte_id in sorted_df_rt_r1.analyte_id.unique():

							analyte_df = sorted_df_rt_r1[sorted_df_rt_r1.analyte_id == analyte_id]

							max_scan_df = analyte_df[analyte_df.intensity == analyte_df.intensity.max()]

							max_peak_id_array = max_scan_df.peak_id.to_numpy()
							max_peak_id = int(max_peak_id_array[0])




							sorted_df_rt1 = analyte_df[analyte_df.peak_id == int(max_peak_id)]


							rt_array_r1.append([sorted_df_rt1[sorted_df_rt1["analyte_id"] == analyte_id]["rt"].to_numpy(),
						                     sorted_df_rt1[sorted_df_rt1["analyte_id"] == analyte_id]["intensity"].to_numpy()])

							analyte_rt_array_r1.append(analyte_id)

					
					analyte_mz_array_r1 = [0 if str(i) == 'nan' else i for i in analyte_mz_array_r1]

					analyte_rt_array_r1 = [0 if str(i) == 'nan' else i for i in analyte_rt_array_r1]



					analyte_colour_mz_array_r1 = [0]*len(analyte_mz_array_r1)
					analyte_colour_rt_array_r1 = [0]*len(analyte_rt_array_r1)

					
					for i in range(0, len(analyte_mz_array_r1)):
						analyte_colour_mz_array_r1[i] = int((analyte_mz_array_r1[i] + 1) % 299)
						


					
					for i in range(0, len(analyte_rt_array_r1)):

						analyte_colour_rt_array_r1[i] = int((analyte_rt_array_r1[i] + 1) % 299)
						
					print('Loading Replicate Tab')


				if chooseAnalyte == "Replicate Analyte ID":



					print('Loading Replicate Tab for Replicate Analyte ID')
					for analyte_id in sorted_df_mz_r1.replicate_analyte_id.unique():

						mz_array_r1.append([sorted_df_mz_r1[sorted_df_mz_r1["replicate_analyte_id"] == analyte_id]["mz"].to_numpy(),
					                     sorted_df_mz_r1[sorted_df_mz_r1["replicate_analyte_id"] == analyte_id]["intensity"].to_numpy(),])

						analyte_mz_array_r1.append(analyte_id)




					print('Loading Replicate Tab for Replicate Analyte ID')

					if Toggle_rt == True:
						print('HERE')

						for analyte_id in sorted_df_rt_r1.replicate_analyte_id.unique():
							rt_array_r1.append([sorted_df_rt_r1[sorted_df_rt_r1["replicate_analyte_id"] == analyte_id]["rt"].to_numpy(),
						                     sorted_df_rt_r1[sorted_df_rt_r1["replicate_analyte_id"] == analyte_id]["intensity"].to_numpy()])
							analyte_rt_array_r1.append(analyte_id)
					else:

						print('HERE')

						for analyte_id in sorted_df_rt_r1.replicate_analyte_id.unique():
							if analyte_id != None:
								print('HERE')

								analyte_df = sorted_df_rt_r1[sorted_df_rt_r1.replicate_analyte_id == analyte_id]
								print('HERE')

								max_scan_df = analyte_df[analyte_df.intensity == analyte_df.intensity.max()]
								print('HERE')
								max_peak_id_array = max_scan_df.peak_id.to_numpy()
								print('HERE')

								max_peak_id = int(max_peak_id_array[0])
								print('HERE')



								sorted_df_rt1 = analyte_df[analyte_df.peak_id == int(max_peak_id)]

			
								rt_array_r1.append([sorted_df_rt1[sorted_df_rt1["replicate_analyte_id"] == analyte_id]["rt"].to_numpy(),
							                     sorted_df_rt1[sorted_df_rt1["replicate_analyte_id"] == analyte_id]["intensity"].to_numpy()])
								analyte_rt_array_r1.append(analyte_id)
							else:
								pass



					analyte_mz_array_r1 = [0 if str(i) == 'nan' else i for i in analyte_mz_array_r1]

					analyte_rt_array_r1 = [0 if str(i) == 'nan' else i for i in analyte_rt_array_r1]

					analyte_mz_array_r1 = [0 if str(i) == 'None' else i for i in analyte_mz_array_r1]

					analyte_rt_array_r1 = [0 if str(i) == 'None' else i for i in analyte_rt_array_r1]

					analyte_colour_mz_array_r1 = [0]*len(analyte_mz_array_r1)
					analyte_colour_rt_array_r1 = [0]*len(analyte_rt_array_r1)

					i = 0
					while i < len(analyte_mz_array_r1):
						analyte_colour_mz_array_r1[i] = int((analyte_mz_array_r1[i] + 1) % 299)

						i+=1


					i = 0
					while i < len(analyte_rt_array_r1):

						analyte_colour_rt_array_r1[i] = int((analyte_rt_array_r1[i] + 1) % 299)
						i+=1




					print('Loading Replicate Tab')




				self.rt_r1_plot.clear()
				self.mz_r1_plot.clear()




				if Toggle_rt == False:



					if comboAnalyte != "Show All":
						print('Loading Replicate Tab')
						alpha=1


						colour=0
						for analyte in mz_array_r1:



							if analyte_mz_array_r1[colour] == float(comboAnalyte):


								self.test = self.mz_r1_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r1[colour],  values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_mz_array_r1[colour],values=5,alpha=255),symbol='o', symbolSize=3)
								self.test.sigPointsClicked.connect(self.clicked)
							else:
								pass


							self.test = self.mz_r1_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r1[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_mz_array_r1[colour],values=5,alpha=alpha),symbol='o', symbolSize=3)
							self.test.sigPointsClicked.connect(self.clicked)
							colour+=1

						colour=0
						for analyte in rt_array_r1:


							if analyte_rt_array_r1[colour] == float(comboAnalyte):


								self.rt_r1_plot.plot( x=analyte[0], y=analyte[1], fillLevel=-0.3,  pen=pg.intColor(analyte_colour_rt_array_r1[colour],  values=5, alpha=255) ,brush=pg.intColor(analyte_colour_rt_array_r1[colour], alpha=50))
								

							else:
								pass

							self.rt_r1_plot.plot( x=analyte[0], y=analyte[1], fillLevel=-0.3, pen=pg.intColor(analyte_colour_rt_array_r1[colour],  values=5, alpha=alpha),brush=pg.intColor(analyte_colour_rt_array_r1[colour], alpha=alpha))
							colour+=1



					else:
						alpha=255
						colour=0
						for analyte in mz_array_r1:


							self.test = self.mz_r1_plot.plot( title="test", x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r1[colour],  values=5), symbolBrush=pg.intColor( analyte_colour_mz_array_r1[colour],values=5),symbol='o', symbolSize=3)
							#test.setAlpha(0, False)
							self.test.sigPointsClicked.connect(self.clicked)
							colour+=1
						colour=0
						for analyte in rt_array_r1:





							self.rt_r1_plot.plot( x=analyte[0], y=analyte[1], fillLevel=-0.3, pen=pg.intColor(analyte_colour_rt_array_r1[colour],  values=5),brush=pg.intColor(analyte_colour_rt_array_r1[colour], alpha=50))
							colour+=1
						print('Loading Replicate Tab')
				else:

					if comboAnalyte != "Show All":
						print('Loading Replicate Tab')
						alpha=1


						colour=0
						for analyte in mz_array_r1:



							if analyte_mz_array_r1[colour] == float(comboAnalyte):


								self.test = self.mz_r1_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r1[colour],  values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_mz_array_r1[colour],values=5,alpha=255),symbol='o', symbolSize=3)
								self.test.sigPointsClicked.connect(self.clicked)
							else:
								pass


							self.test = self.mz_r1_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r1[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_mz_array_r1[colour],values=5,alpha=alpha),symbol='o', symbolSize=3)
							self.test.sigPointsClicked.connect(self.clicked)
							colour+=1

						colour=0
						for analyte in rt_array_r1:


							if analyte_rt_array_r1[colour] == float(comboAnalyte):


								self.rt_r1_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r1[colour],  values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_rt_array_r1[colour], values=5, alpha=255),symbol='o', symbolSize=3)
								

							else:
								pass

							self.rt_r1_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r1[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_rt_array_r1[colour], values=5, alpha=alpha),symbol='o', symbolSize=3)
							colour+=1



					else:
						alpha=255
						colour=0
						for analyte in mz_array_r1:


							self.test = self.mz_r1_plot.plot( title="test", x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r1[colour],  values=5), symbolBrush=pg.intColor( analyte_colour_mz_array_r1[colour],values=5),symbol='o', symbolSize=3)
							#test.setAlpha(0, False)
							self.test.sigPointsClicked.connect(self.clicked)
							colour+=1
						colour=0
						for analyte in rt_array_r1:





							self.rt_r1_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r1[colour],  values=5), symbolBrush=pg.intColor( analyte_colour_rt_array_r1[colour], values=5),symbol='o', symbolSize=3)
							colour+=1
						print('Loading Replicate Tab')

				############### Fetch data for rt plot and plot sample from comboBox ###############




				mz_array_r2 = []
				analyte_mz_array_r2 = []
				rt_array_r2 = []
				analyte_rt_array_r2 = []

				if chooseAnalyte == "Sample Analyte ID":
					print('Loading Replicate Tab')
					for analyte_id in sorted_df_mz_r2.analyte_id.unique():

						mz_array_r2.append([sorted_df_mz_r2[sorted_df_mz_r2["analyte_id"] == analyte_id]["mz"].to_numpy(),
					                     sorted_df_mz_r2[sorted_df_mz_r2["analyte_id"] == analyte_id]["intensity"].to_numpy(),])

						analyte_mz_array_r2.append(analyte_id)


						



					if Toggle_rt == True:
					

						for analyte_id in sorted_df_rt_r1.analyte_id.unique():
							rt_array_r2.append([sorted_df_rt_r2[sorted_df_rt_r2["analyte_id"] == analyte_id]["rt"].to_numpy(),
						                     sorted_df_rt_r2[sorted_df_rt_r2["analyte_id"] == analyte_id]["intensity"].to_numpy()])
							analyte_rt_array_r2.append(analyte_id)
					else:



						for analyte_id in sorted_df_rt_r2.analyte_id.unique():

							analyte_df = sorted_df_rt_r2[sorted_df_rt_r2.analyte_id == analyte_id]

							max_scan_df = analyte_df[analyte_df.intensity == analyte_df.intensity.max()]

							max_peak_id_array = max_scan_df.peak_id.to_numpy()
							max_peak_id = int(max_peak_id_array[0])


		


							sorted_df_rt1 = analyte_df[analyte_df.peak_id == int(max_peak_id)]


							rt_array_r2.append([sorted_df_rt1[sorted_df_rt1["analyte_id"] == analyte_id]["rt"].to_numpy(),
						                     sorted_df_rt1[sorted_df_rt1["analyte_id"] == analyte_id]["intensity"].to_numpy()])
							analyte_rt_array_r2.append(analyte_id)


					analyte_mz_array_r2 = [0 if str(i) == 'nan' else i for i in analyte_mz_array_r2]

					analyte_rt_array_r2 = [0 if str(i) == 'nan' else i for i in analyte_rt_array_r2]



					analyte_colour_mz_array_r2 = [0]*len(analyte_mz_array_r2)
					analyte_colour_rt_array_r2 = [0]*len(analyte_rt_array_r2)
					print('Loading Replicate Tab')
					i = 0
					while i < len(analyte_mz_array_r2):
						analyte_colour_mz_array_r2[i] = int((analyte_mz_array_r2[i] + 1) % 299)
						i+=1


					i = 0
					while i < len(analyte_rt_array_r2):

						analyte_colour_rt_array_r2[i] = int((analyte_rt_array_r2[i] + 1) % 299)
						i+=1


				if chooseAnalyte == "Replicate Analyte ID":
					print('Loading Replicate Tab')
					for analyte_id in sorted_df_mz_r2.replicate_analyte_id.unique():

						mz_array_r2.append([sorted_df_mz_r2[sorted_df_mz_r2["replicate_analyte_id"] == analyte_id]["mz"].to_numpy(),
					                     sorted_df_mz_r2[sorted_df_mz_r2["replicate_analyte_id"] == analyte_id]["intensity"].to_numpy(),])

						analyte_mz_array_r2.append(analyte_id)


						



					if Toggle_rt == True:
					

						for analyte_id in sorted_df_rt_r1.replicate_analyte_id.unique():
							rt_array_r2.append([sorted_df_rt_r2[sorted_df_rt_r2["replicate_analyte_id"] == analyte_id]["rt"].to_numpy(),
						                     sorted_df_rt_r2[sorted_df_rt_r2["replicate_analyte_id"] == analyte_id]["intensity"].to_numpy()])
							analyte_rt_array_r2.append(analyte_id)
					else:



						for analyte_id in sorted_df_rt_r2.replicate_analyte_id.unique():
							if analyte_id != None:
								analyte_df = sorted_df_rt_r2[sorted_df_rt_r2.replicate_analyte_id == analyte_id]

								max_scan_df = analyte_df[analyte_df.intensity == analyte_df.intensity.max()]

								max_peak_id_array = max_scan_df.peak_id.to_numpy()
								max_peak_id = int(max_peak_id_array[0])




								sorted_df_rt1 = analyte_df[analyte_df.peak_id == int(max_peak_id)]

			
								rt_array_r2.append([sorted_df_rt1[sorted_df_rt1["replicate_analyte_id"] == analyte_id]["rt"].to_numpy(),
							                     sorted_df_rt1[sorted_df_rt1["replicate_analyte_id"] == analyte_id]["intensity"].to_numpy()])
								analyte_rt_array_r2.append(analyte_id)
							else:
								pass


					analyte_mz_array_r2 = [0 if str(i) == 'nan' else i for i in analyte_mz_array_r2]

					analyte_rt_array_r2 = [0 if str(i) == 'nan' else i for i in analyte_rt_array_r2]

					analyte_mz_array_r2 = [0 if str(i) == 'None' else i for i in analyte_mz_array_r2]

					analyte_rt_array_r2 = [0 if str(i) == 'None' else i for i in analyte_rt_array_r2]

					analyte_colour_mz_array_r2 = [0]*len(analyte_mz_array_r2)
					analyte_colour_rt_array_r2 = [0]*len(analyte_rt_array_r2)

					i = 0
					while i < len(analyte_mz_array_r2):
						analyte_colour_mz_array_r2[i] = int((analyte_mz_array_r2[i] + 1) % 299)
						i+=1


					i = 0
					while i < len(analyte_rt_array_r2):

						analyte_colour_rt_array_r2[i] = int((analyte_rt_array_r2[i] + 1) % 299)
						i+=1





				if chooseAnalyte == "Experiment Analyte ID":
					print('Loading Replicate Tab')
					for analyte_id in sorted_df_mz_r2.experiment_analyte_id.unique():

						mz_array_r2.append([sorted_df_mz_r2[sorted_df_mz_r2["experiment_analyte_id"] == analyte_id]["mz"].to_numpy(),
					                     sorted_df_mz_r2[sorted_df_mz_r2["experiment_analyte_id"] == analyte_id]["intensity"].to_numpy(),])

						analyte_mz_array_r2.append(analyte_id)


						



					if Toggle_rt == True:
					

						for analyte_id in sorted_df_rt_r1.experiment_analyte_id.unique():
							rt_array_r2.append([sorted_df_rt_r2[sorted_df_rt_r2["experiment_analyte_id"] == analyte_id]["rt"].to_numpy(),
						                     sorted_df_rt_r2[sorted_df_rt_r2["experiment_analyte_id"] == analyte_id]["intensity"].to_numpy()])
							analyte_rt_array_r2.append(analyte_id)
					else:



						for analyte_id in sorted_df_rt_r2.experiment_analyte_id.unique():
							if analyte_id != None:
								analyte_df = sorted_df_rt_r2[sorted_df_rt_r2.experiment_analyte_id == analyte_id]

								max_scan_df = analyte_df[analyte_df.intensity == analyte_df.intensity.max()]

								max_peak_id_array = max_scan_df.peak_id.to_numpy()
								max_peak_id = int(max_peak_id_array[0])




								sorted_df_rt1 = analyte_df[analyte_df.peak_id == int(max_peak_id)]

			
								rt_array_r2.append([sorted_df_rt1[sorted_df_rt1["experiment_analyte_id"] == analyte_id]["rt"].to_numpy(),
							                     sorted_df_rt1[sorted_df_rt1["experiment_analyte_id"] == analyte_id]["intensity"].to_numpy()])
								analyte_rt_array_r2.append(analyte_id)
							else:
								pass


					analyte_mz_array_r2 = [0 if str(i) == 'nan' else i for i in analyte_mz_array_r2]

					analyte_rt_array_r2 = [0 if str(i) == 'nan' else i for i in analyte_rt_array_r2]

					analyte_mz_array_r2 = [0 if str(i) == 'None' else i for i in analyte_mz_array_r2]

					analyte_rt_array_r2 = [0 if str(i) == 'None' else i for i in analyte_rt_array_r2]

					analyte_colour_mz_array_r2 = [0]*len(analyte_mz_array_r2)
					analyte_colour_rt_array_r2 = [0]*len(analyte_rt_array_r2)

					i = 0
					while i < len(analyte_mz_array_r2):
						analyte_colour_mz_array_r2[i] = int((analyte_mz_array_r2[i] + 1) % 299)
						i+=1


					i = 0
					while i < len(analyte_rt_array_r2):

						analyte_colour_rt_array_r2[i] = int((analyte_rt_array_r2[i] + 1) % 299)
						i+=1











				self.rt_r2_plot.clear()
				self.mz_r2_plot.clear()
				print('Loading Replicate Tab')


				if Toggle_rt == False:



					if comboAnalyte != "Show All":
						print('Loading Replicate Tab')
						alpha=1


						colour=0
						for analyte in mz_array_r2:


							#self.mz_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_array[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_array[colour],values=5,alpha=alpha),symbol='o', symbolSize=3)
							#colour+=1

							if analyte_mz_array_r2[colour] == float(comboAnalyte):


								self.test = self.mz_r2_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r2[colour],  values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_mz_array_r2[colour],values=5,alpha=255),symbol='o', symbolSize=3)
								self.test.sigPointsClicked.connect(self.clicked)
							else:
								pass


							self.test = self.mz_r2_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r2[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_mz_array_r2[colour],values=5,alpha=alpha),symbol='o', symbolSize=3)
							self.test.sigPointsClicked.connect(self.clicked)
							colour+=1

						colour=0
						for analyte in rt_array_r2:


							if analyte_rt_array_r2[colour] == float(comboAnalyte):


								self.rt_r2_plot.plot( x=analyte[0], y=analyte[1], fillLevel=-0.3, pen=pg.intColor(analyte_colour_rt_array_r2[colour],  values=5, alpha=255),brush=pg.intColor(analyte_colour_rt_array_r2[colour], alpha=50))
								

							else:
								pass

							self.rt_r2_plot.plot( x=analyte[0], y=analyte[1], fillLevel=-0.3, pen=pg.intColor(analyte_colour_rt_array_r2[colour],  values=5, alpha=alpha),brush=pg.intColor(analyte_colour_rt_array_r2[colour], alpha=alpha))
							colour+=1


						print('Loading Replicate Tab')
					else:
						alpha=255
						colour=0
						for analyte in mz_array_r2:


							self.test = self.mz_r2_plot.plot( title="test", x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r2[colour],  values=5), symbolBrush=pg.intColor( analyte_colour_mz_array_r2[colour],values=5),symbol='o', symbolSize=3)
							#test.setAlpha(0, False)
							self.test.sigPointsClicked.connect(self.clicked)
							colour+=1
						colour=0
						for analyte in rt_array_r2:





							self.rt_r2_plot.plot( x=analyte[0], y=analyte[1], fillLevel=-0.3, pen=pg.intColor(analyte_colour_rt_array_r2[colour],  values=5),brush=pg.intColor(analyte_colour_rt_array_r2[colour], alpha=50))
							colour+=1
				else:




					if comboAnalyte != "Show All":
							print('Loading Replicate Tab')
							alpha=1


							colour=0
							for analyte in mz_array_r2:


								#self.mz_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_array[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_array[colour],values=5,alpha=alpha),symbol='o', symbolSize=3)
								#colour+=1

								if analyte_mz_array_r2[colour] == float(comboAnalyte):


									self.test = self.mz_r2_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r2[colour],  values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_mz_array_r2[colour],values=5,alpha=255),symbol='o', symbolSize=3)
									self.test.sigPointsClicked.connect(self.clicked)
								else:
									pass


								self.test = self.mz_r2_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r2[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_mz_array_r2[colour],values=5,alpha=alpha),symbol='o', symbolSize=3)
								self.test.sigPointsClicked.connect(self.clicked)
								colour+=1

							colour=0
							for analyte in rt_array_r2:


								if analyte_rt_array_r2[colour] == float(comboAnalyte):


									self.rt_r2_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r2[colour],  values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_rt_array_r2[colour], values=5, alpha=255),symbol='o', symbolSize=3)
									

								else:
									pass

								self.rt_r2_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r2[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_rt_array_r2[colour], values=5, alpha=alpha),symbol='o', symbolSize=3)
								colour+=1


							print('Loading Replicate Tab')
					else:



						alpha=255
						colour=0
						for analyte in mz_array_r2:


							self.test = self.mz_r2_plot.plot( title="test", x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r2[colour],  values=5), symbolBrush=pg.intColor( analyte_colour_mz_array_r2[colour],values=5),symbol='o', symbolSize=3)
							#test.setAlpha(0, False)
							self.test.sigPointsClicked.connect(self.clicked)
							colour+=1
						colour=0
						for analyte in rt_array_r2:





							self.rt_r2_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r2[colour],  values=5), symbolBrush=pg.intColor( analyte_colour_rt_array_r2[colour], values=5),symbol='o', symbolSize=3)
							colour+=1



				mz_array_r3 = []
				analyte_mz_array_r3 = []
				rt_array_r3 = []
				analyte_rt_array_r3 = []
				print('Loading Replicate Tab')





				if chooseAnalyte == "Sample Analyte ID":
					print('Loading Replicate Tab')
					for analyte_id in sorted_df_mz_r3.analyte_id.unique():

						mz_array_r3.append([sorted_df_mz_r3[sorted_df_mz_r3["analyte_id"] == analyte_id]["mz"].to_numpy(),
					                     sorted_df_mz_r3[sorted_df_mz_r3["analyte_id"] == analyte_id]["intensity"].to_numpy(),])

						analyte_mz_array_r3.append(analyte_id)


						



					if Toggle_rt == True:
					

						for analyte_id in sorted_df_rt_r3.analyte_id.unique():
							rt_array_r3.append([sorted_df_rt_r3[sorted_df_rt_r3["analyte_id"] == analyte_id]["rt"].to_numpy(),
						                     sorted_df_rt_r3[sorted_df_rt_r3["analyte_id"] == analyte_id]["intensity"].to_numpy()])
							analyte_rt_array_r3.append(analyte_id)
					else:



						for analyte_id in sorted_df_rt_r3.analyte_id.unique():

							analyte_df = sorted_df_rt_r3[sorted_df_rt_r3.analyte_id == analyte_id]

							max_scan_df = analyte_df[analyte_df.intensity == analyte_df.intensity.max()]

							max_peak_id_array = max_scan_df.peak_id.to_numpy()
							max_peak_id = int(max_peak_id_array[0])

		


							sorted_df_rt1 = analyte_df[analyte_df.peak_id == int(max_peak_id)]

		
							rt_array_r3.append([sorted_df_rt1[sorted_df_rt1["analyte_id"] == analyte_id]["rt"].to_numpy(),
						                     sorted_df_rt1[sorted_df_rt1["analyte_id"] == analyte_id]["intensity"].to_numpy()])
							analyte_rt_array_r3.append(analyte_id)


					analyte_mz_array_r3 = [0 if str(i) == 'nan' else i for i in analyte_mz_array_r3]

					analyte_rt_array_r3 = [0 if str(i) == 'nan' else i for i in analyte_rt_array_r3]



					analyte_colour_mz_array_r3 = [0]*len(analyte_mz_array_r3)
					analyte_colour_rt_array_r3 = [0]*len(analyte_rt_array_r3)

					i = 0
					while i < len(analyte_mz_array_r3):
						analyte_colour_mz_array_r3[i] = int((analyte_mz_array_r3[i] + 1) % 299)
						i+=1


					i = 0
					while i < len(analyte_rt_array_r3):

						analyte_colour_rt_array_r3[i] = int((analyte_rt_array_r3[i] + 1) % 299)
						i+=1


				if chooseAnalyte == "Replicate Analyte ID":



					for analyte_id in sorted_df_mz_r3.replicate_analyte_id.unique():

						mz_array_r3.append([sorted_df_mz_r3[sorted_df_mz_r3["replicate_analyte_id"] == analyte_id]["mz"].to_numpy(),
					                     sorted_df_mz_r3[sorted_df_mz_r3["replicate_analyte_id"] == analyte_id]["intensity"].to_numpy(),])

						analyte_mz_array_r3.append(analyte_id)


						

					print('Loading Replicate Tab')

					if Toggle_rt == True:
					

						for analyte_id in sorted_df_rt_r3.replicate_analyte_id.unique():
							rt_array_r3.append([sorted_df_rt_r3[sorted_df_rt_r3["replicate_analyte_id"] == analyte_id]["rt"].to_numpy(),
						                     sorted_df_rt_r3[sorted_df_rt_r3["replicate_analyte_id"] == analyte_id]["intensity"].to_numpy()])
							analyte_rt_array_r3.append(analyte_id)
					else:



						for analyte_id in sorted_df_rt_r3.replicate_analyte_id.unique():
							if analyte_id != None:
								analyte_df = sorted_df_rt_r3[sorted_df_rt_r3.replicate_analyte_id == analyte_id]

								max_scan_df = analyte_df[analyte_df.intensity == analyte_df.intensity.max()]

								max_peak_id_array = max_scan_df.peak_id.to_numpy()
								max_peak_id = int(max_peak_id_array[0])

			


								sorted_df_rt1 = analyte_df[analyte_df.peak_id == int(max_peak_id)]

			
								rt_array_r3.append([sorted_df_rt1[sorted_df_rt1["replicate_analyte_id"] == analyte_id]["rt"].to_numpy(),
							                     sorted_df_rt1[sorted_df_rt1["replicate_analyte_id"] == analyte_id]["intensity"].to_numpy()])
								analyte_rt_array_r3.append(analyte_id)
							else:
								pass


					analyte_mz_array_r3 = [0 if str(i) == 'nan' else i for i in analyte_mz_array_r3]

					analyte_rt_array_r3 = [0 if str(i) == 'nan' else i for i in analyte_rt_array_r3]

					analyte_mz_array_r3 = [0 if str(i) == 'None' else i for i in analyte_mz_array_r3]

					analyte_rt_array_r3 = [0 if str(i) == 'None' else i for i in analyte_rt_array_r3]

					analyte_colour_mz_array_r3 = [0]*len(analyte_mz_array_r3)
					analyte_colour_rt_array_r3 = [0]*len(analyte_rt_array_r3)

					i = 0
					while i < len(analyte_mz_array_r3):
						analyte_colour_mz_array_r3[i] = int((analyte_mz_array_r3[i] + 1) % 299)
						i+=1


					i = 0
					while i < len(analyte_rt_array_r3):

						analyte_colour_rt_array_r3[i] = int((analyte_rt_array_r3[i] + 1) % 299)
						i+=1


				if chooseAnalyte == "Experiment Analyte ID":


					for analyte_id in sorted_df_mz_r3.experiment_analyte_id.unique():

						mz_array_r3.append([sorted_df_mz_r3[sorted_df_mz_r3["experiment_analyte_id"] == analyte_id]["mz"].to_numpy(),
					                     sorted_df_mz_r3[sorted_df_mz_r3["experiment_analyte_id"] == analyte_id]["intensity"].to_numpy(),])

						analyte_mz_array_r3.append(analyte_id)


						

					print('Loading Replicate Tab')

					if Toggle_rt == True:
					

						for analyte_id in sorted_df_rt_r3.experiment_analyte_id.unique():
							rt_array_r3.append([sorted_df_rt_r3[sorted_df_rt_r3["experiment_analyte_id"] == analyte_id]["rt"].to_numpy(),
						                     sorted_df_rt_r3[sorted_df_rt_r3["experiment_analyte_id"] == analyte_id]["intensity"].to_numpy()])
							analyte_rt_array_r3.append(analyte_id)
					else:



						for analyte_id in sorted_df_rt_r3.experiment_analyte_id.unique():
							if analyte_id != None:
								analyte_df = sorted_df_rt_r3[sorted_df_rt_r3.experiment_analyte_id == analyte_id]

								max_scan_df = analyte_df[analyte_df.intensity == analyte_df.intensity.max()]

								max_peak_id_array = max_scan_df.peak_id.to_numpy()
								max_peak_id = int(max_peak_id_array[0])

			


								sorted_df_rt1 = analyte_df[analyte_df.peak_id == int(max_peak_id)]

			
								rt_array_r3.append([sorted_df_rt1[sorted_df_rt1["experiment_analyte_id"] == analyte_id]["rt"].to_numpy(),
							                     sorted_df_rt1[sorted_df_rt1["experiment_analyte_id"] == analyte_id]["intensity"].to_numpy()])
								analyte_rt_array_r3.append(analyte_id)
							else:
								pass


					analyte_mz_array_r3 = [0 if str(i) == 'nan' else i for i in analyte_mz_array_r3]

					analyte_rt_array_r3 = [0 if str(i) == 'nan' else i for i in analyte_rt_array_r3]

					analyte_mz_array_r3 = [0 if str(i) == 'None' else i for i in analyte_mz_array_r3]

					analyte_rt_array_r3 = [0 if str(i) == 'None' else i for i in analyte_rt_array_r3]

					analyte_colour_mz_array_r3 = [0]*len(analyte_mz_array_r3)
					analyte_colour_rt_array_r3 = [0]*len(analyte_rt_array_r3)

					i = 0
					while i < len(analyte_mz_array_r3):
						analyte_colour_mz_array_r3[i] = int((analyte_mz_array_r3[i] + 1) % 299)
						i+=1


					i = 0
					while i < len(analyte_rt_array_r3):

						analyte_colour_rt_array_r3[i] = int((analyte_rt_array_r3[i] + 1) % 299)
						i+=1






				self.rt_r3_plot.clear()
				self.mz_r3_plot.clear()
				print('Loading Replicate Tab')



				if Toggle_rt == False:
					if comboAnalyte != "Show All":
						print('Loading Replicate Tab')
						alpha=1


						colour=0
						for analyte in mz_array_r3:


							#self.mz_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_array[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_array[colour],values=5,alpha=alpha),symbol='o', symbolSize=3)
							#colour+=1

							if analyte_mz_array_r3[colour] == float(comboAnalyte):


								self.test = self.mz_r3_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r3[colour],  values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_mz_array_r3[colour],values=5,alpha=255),symbol='o', symbolSize=3)
								self.test.sigPointsClicked.connect(self.clicked)
							else:
								pass


							self.test = self.mz_r3_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r3[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_mz_array_r3[colour],values=5,alpha=alpha),symbol='o', symbolSize=3)
							self.test.sigPointsClicked.connect(self.clicked)
							colour+=1

						colour=0
						for analyte in rt_array_r3:


							if analyte_rt_array_r3[colour] == float(comboAnalyte):


								self.rt_r3_plot.plot( x=analyte[0], y=analyte[1], fillLevel=-0.3, pen=pg.intColor(analyte_colour_rt_array_r3[colour],  values=5, alpha=255),brush=pg.intColor(analyte_colour_rt_array_r3[colour], alpha=50))
								

							else:
								pass

							self.rt_r3_plot.plot( x=analyte[0], y=analyte[1], fillLevel=-0.3, pen=pg.intColor(analyte_colour_rt_array_r3[colour],  values=5, alpha=alpha),brush=pg.intColor(analyte_colour_rt_array_r3[colour], alpha=alpha))
							colour+=1



					else:
						alpha=255
						colour=0
						for analyte in mz_array_r3:


							self.test = self.mz_r3_plot.plot( title="test", x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r3[colour],  values=5), symbolBrush=pg.intColor( analyte_colour_mz_array_r3[colour],values=5),symbol='o', symbolSize=3)
							#test.setAlpha(0, False)
							self.test.sigPointsClicked.connect(self.clicked)
							colour+=1
						colour=0
						for analyte in rt_array_r3:





							self.rt_r3_plot.plot( x=analyte[0], y=analyte[1], fillLevel=-0.3, pen=pg.intColor(analyte_colour_rt_array_r3[colour],  values=5),brush=pg.intColor(analyte_colour_rt_array_r3[colour], alpha=50))
							colour+=1
				else:
					if comboAnalyte != "Show All":
							print('Loading Replicate Tab')
							alpha=1


							colour=0
							for analyte in mz_array_r3:


								#self.mz_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_array[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_array[colour],values=5,alpha=alpha),symbol='o', symbolSize=3)
								#colour+=1

								if analyte_mz_array_r3[colour] == float(comboAnalyte):


									self.test = self.mz_r3_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r3[colour],  values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_mz_array_r3[colour],values=5,alpha=255),symbol='o', symbolSize=3)
									self.test.sigPointsClicked.connect(self.clicked)
								else:
									pass


								self.test = self.mz_r3_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r3[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_mz_array_r3[colour],values=5,alpha=alpha),symbol='o', symbolSize=3)
								self.test.sigPointsClicked.connect(self.clicked)
								colour+=1

							colour=0
							for analyte in rt_array_r3:


								if analyte_rt_array_r3[colour] == float(comboAnalyte):


									self.rt_r3_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r3[colour],  values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_rt_array_r3[colour], values=5, alpha=255),symbol='o', symbolSize=3)
									

								else:
									pass

								self.rt_r3_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r3[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_rt_array_r3[colour], values=5, alpha=alpha),symbol='o', symbolSize=3)
								colour+=1



					else:
						alpha=255
						colour=0
						for analyte in mz_array_r3:


							self.test = self.mz_r3_plot.plot( title="test", x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r3[colour],  values=5), symbolBrush=pg.intColor( analyte_colour_mz_array_r3[colour],values=5),symbol='o', symbolSize=3)
							#test.setAlpha(0, False)
							self.test.sigPointsClicked.connect(self.clicked)
							colour+=1
						colour=0
						for analyte in rt_array_r3:





							self.rt_r3_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r3[colour],  values=5), symbolBrush=pg.intColor( analyte_colour_rt_array_r3[colour], values=5),symbol='o', symbolSize=3)
							colour+=1



				if chooseAnalyte == "Experiment Analyte ID":


					for analyte_id in sorted_df_mz_r1.experiment_analyte_id.unique():

						mz_array_r1.append([sorted_df_mz_r1[sorted_df_mz_r1["experiment_analyte_id"] == analyte_id]["mz"].to_numpy(),
					                     sorted_df_mz_r1[sorted_df_mz_r1["experiment_analyte_id"] == analyte_id]["intensity"].to_numpy(),])

						analyte_mz_array_r1.append(analyte_id)


						

					print('Loading Replicate Tab')

					if Toggle_rt == True:
					

						for analyte_id in sorted_df_rt_r1.experiment_analyte_id.unique():
							rt_array_r1.append([sorted_df_rt_r1[sorted_df_rt_r1["experiment_analyte_id"] == analyte_id]["rt"].to_numpy(),
						                     sorted_df_rt_r1[sorted_df_rt_r1["experiment_analyte_id"] == analyte_id]["intensity"].to_numpy()])
							analyte_rt_array_r1.append(analyte_id)
					else:



						for analyte_id in sorted_df_rt_r1.experiment_analyte_id.unique():
							if analyte_id != None:
								analyte_df = sorted_df_rt_r1[sorted_df_rt_r1.experiment_analyte_id == analyte_id]

								max_scan_df = analyte_df[analyte_df.intensity == analyte_df.intensity.max()]

								max_peak_id_array = max_scan_df.peak_id.to_numpy()
								max_peak_id = int(max_peak_id_array[0])

			


								sorted_df_rt1 = analyte_df[analyte_df.peak_id == int(max_peak_id)]

			
								rt_array_r1.append([sorted_df_rt1[sorted_df_rt1["experiment_analyte_id"] == analyte_id]["rt"].to_numpy(),
							                     sorted_df_rt1[sorted_df_rt1["experiment_analyte_id"] == analyte_id]["intensity"].to_numpy()])
								analyte_rt_array_r1.append(analyte_id)
							else:
								pass


					analyte_mz_array_r1 = [0 if str(i) == 'nan' else i for i in analyte_mz_array_r1]

					analyte_rt_array_r1 = [0 if str(i) == 'nan' else i for i in analyte_rt_array_r1]

					analyte_mz_array_r1 = [0 if str(i) == 'None' else i for i in analyte_mz_array_r1]

					analyte_rt_array_r1 = [0 if str(i) == 'None' else i for i in analyte_rt_array_r1]

					analyte_colour_mz_array_r1 = [0]*len(analyte_mz_array_r1)
					analyte_colour_rt_array_r1 = [0]*len(analyte_rt_array_r1)

					i = 0
					while i < len(analyte_mz_array_r1):
						analyte_colour_mz_array_r1[i] = int((analyte_mz_array_r1[i] + 1) % 299)
						i+=1


					i = 0
					while i < len(analyte_rt_array_r1):

						analyte_colour_rt_array_r1[i] = int((analyte_rt_array_r1[i] + 1) % 299)
						i+=1






				self.rt_r1_plot.clear()
				self.mz_r1_plot.clear()
				print('Loading Replicate Tab')



				if Toggle_rt == False:
					if comboAnalyte != "Show All":
						print('Loading Replicate Tab')
						alpha=1


						colour=0
						for analyte in mz_array_r1:


							#self.mz_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_array[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_array[colour],values=5,alpha=alpha),symbol='o', symbolSize=3)
							#colour+=1

							if analyte_mz_array_r1[colour] == float(comboAnalyte):


								self.test = self.mz_r1_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r1[colour],  values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_mz_array_r1[colour],values=5,alpha=255),symbol='o', symbolSize=3)
								self.test.sigPointsClicked.connect(self.clicked)
							else:
								pass


							self.test = self.mz_r1_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r1[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_mz_array_r1[colour],values=5,alpha=alpha),symbol='o', symbolSize=3)
							self.test.sigPointsClicked.connect(self.clicked)
							colour+=1

						colour=0
						for analyte in rt_array_r1:


							if analyte_rt_array_r1[colour] == float(comboAnalyte):


								self.rt_r1_plot.plot( x=analyte[0], y=analyte[1], fillLevel=-0.3, pen=pg.intColor(analyte_colour_rt_array_r1[colour],  values=5, alpha=255),brush=pg.intColor(analyte_colour_rt_array_r1[colour], alpha=50))
								

							else:
								pass

							self.rt_r1_plot.plot( x=analyte[0], y=analyte[1], fillLevel=-0.3, pen=pg.intColor(analyte_colour_rt_array_r1[colour],  values=5, alpha=alpha),brush=pg.intColor(analyte_colour_rt_array_r1[colour], alpha=alpha))
							colour+=1



					else:
						alpha=255
						colour=0
						for analyte in mz_array_r1:


							self.test = self.mz_r1_plot.plot( title="test", x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r1[colour],  values=5), symbolBrush=pg.intColor( analyte_colour_mz_array_r1[colour],values=5),symbol='o', symbolSize=3)
							#test.setAlpha(0, False)
							self.test.sigPointsClicked.connect(self.clicked)
							colour+=1
						colour=0
						for analyte in rt_array_r1:





							self.rt_r1_plot.plot( x=analyte[0], y=analyte[1], fillLevel=-0.3, pen=pg.intColor(analyte_colour_rt_array_r1[colour],  values=5),brush=pg.intColor(analyte_colour_rt_array_r1[colour], alpha=50))
							colour+=1
				else:
					if comboAnalyte != "Show All":
							print('Loading Replicate Tab')
							alpha=1


							colour=0
							for analyte in mz_array_r1:


								#self.mz_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_array[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_array[colour],values=5,alpha=alpha),symbol='o', symbolSize=3)
								#colour+=1

								if analyte_mz_array_r1[colour] == float(comboAnalyte):


									self.test = self.mz_r1_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r1[colour],  values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_mz_array_r1[colour],values=5,alpha=255),symbol='o', symbolSize=3)
									self.test.sigPointsClicked.connect(self.clicked)
								else:
									pass


								self.test = self.mz_r1_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r1[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_mz_array_r1[colour],values=5,alpha=alpha),symbol='o', symbolSize=3)
								self.test.sigPointsClicked.connect(self.clicked)
								colour+=1

							colour=0
							for analyte in rt_array_r1:


								if analyte_rt_array_r1[colour] == float(comboAnalyte):


									self.rt_r1_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r1[colour],  values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_rt_array_r1[colour], values=5, alpha=255),symbol='o', symbolSize=3)
									

								else:
									pass

								self.rt_r1_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r1[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_rt_array_r1[colour], values=5, alpha=alpha),symbol='o', symbolSize=3)
								colour+=1



					else:
						alpha=255
						colour=0
						for analyte in mz_array_r1:


							self.test = self.mz_r1_plot.plot( title="test", x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r1[colour],  values=5), symbolBrush=pg.intColor( analyte_colour_mz_array_r1[colour],values=5),symbol='o', symbolSize=3)
							#test.setAlpha(0, False)
							self.test.sigPointsClicked.connect(self.clicked)
							colour+=1
						colour=0
						for analyte in rt_array_r1:





							self.rt_r1_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r1[colour],  values=5), symbolBrush=pg.intColor( analyte_colour_rt_array_r1[colour], values=5),symbol='o', symbolSize=3)
							colour+=1
				#analyte_array[0] = 0
				#print(analyte_array)

				self.mz_r1_plot.getAxis('left').setStyle(showValues=False)
				self.mz_r2_plot.getAxis('left').setStyle(showValues=False)
				self.mz_r3_plot.getAxis('left').setStyle(showValues=False)
				self.rt_r1_plot.getAxis('left').setStyle(showValues=False)
				self.rt_r2_plot.getAxis('left').setStyle(showValues=False)
				self.rt_r3_plot.getAxis('left').setStyle(showValues=False)




				self.mz_r1_plot.enableAutoRange(axis=None)
				self.mz_r2_plot.enableAutoRange(axis=None)
				self.mz_r3_plot.enableAutoRange(axis=None)
				self.rt_r1_plot.enableAutoRange(axis=None)
				self.rt_r2_plot.enableAutoRange(axis=None)
				self.rt_r3_plot.enableAutoRange(axis=None)
				state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), self.rtMinBox.text(), self.rtMaxBox.text(),self.Toggle_rt.isChecked())

				self.mz_r1_plot.setXRange(float(state.return_massMin_state()),float(state.return_massMax_state()))
				self.mz_r2_plot.setXRange(float(state.return_massMin_state()),float(state.return_massMax_state()))
				self.mz_r3_plot.setXRange(float(state.return_massMin_state()),float(state.return_massMax_state()))
				self.rt_r1_plot.setXRange(float(state.return_rtMin_state()),float(state.return_rtMax_state()))
				self.rt_r2_plot.setXRange(float(state.return_rtMin_state()),float(state.return_rtMax_state()))
				self.rt_r3_plot.setXRange(float(state.return_rtMin_state()),float(state.return_rtMax_state()))
			else:
				pass
			self.c.append('Replicate Tab Loaded Succesfully')
			print('PLOT REPLICATE END')
		return

	def Plot_Experiment(self):
		if open_state ==2:
			print('PLOT EXPERIMENT BEGIN')
			# self.ex_bOn_plot.setYRange(min=-20,max=1)
			# self.ex_bOn_plot.setXRange(min=100,max=1200)
			self.ex_bOn_plot.disableAutoRange(axis=None)
			tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
			Sample_State = tab_state.return_sample_state()
			Replicate_State = tab_state.return_replicate_state()
			Experiment_State = tab_state.return_experiment_state()



			if Experiment_State == False:
				comboAnalyte=self.comboAnalyte.currentText()

				

				samples, experiment_df = self.Experiment_Plot_DF()

				self.ex_bOn_plot.clear()
				
				for i in range(0, len(samples)):
					analyte_array = []
					experiment_sample_df = experiment_df[experiment_df.sample_name == samples[i]]


					experiment_array_BOn = []


					for experiment_analyte_id in experiment_sample_df.experiment_analyte_id.unique():
						experiment_array_BOn.append([experiment_sample_df[experiment_sample_df["experiment_analyte_id"] == experiment_analyte_id]["experiment_analyte_max_mass"].to_numpy(), 
										experiment_sample_df[experiment_sample_df["experiment_analyte_id"] == experiment_analyte_id]["sample_analyte_max_intensity"].to_numpy()])
						analyte_array.append(experiment_analyte_id)

					
					analyte_array = [0 if str(z) == 'nan' else z for z in analyte_array]



					analyte_colour_array = [0]*len(analyte_array)


					
					for a in range(0, len(analyte_array)):

						analyte_colour_array[a] = int((analyte_array[a] + 1) % 299)
						

					if comboAnalyte == "Show All":
						
						colour = 0
						for analyte in experiment_array_BOn:





							y=[i]*len(analyte[0])

							x=analyte[0]

							max_int = max(analyte[1])


							test = self.ex_bOn_plot.plot(title= "Test", x=x, y=y, pen=None, symbolPen=pg.intColor(analyte_colour_array[colour],  values=5), symbolBrush=pg.intColor(analyte_colour_array[colour],values=5),symbol='o', symbolSize= 10)
							min_int = 50000
							int_x=40
							if max_int < min_int:
								test.setSymbolSize((max_int/230000)*int_x)
							else:
								test.setSymbolSize((max_int/230000)*4)

							colour+=1
					if comboAnalyte != "Show All":
						alpha=30
						colour = 0

						for analyte in experiment_array_BOn:




							max_int = max(analyte[1])
							if analyte_array[colour] == float(comboAnalyte):


								test = self.ex_bOn_plot.plot(title= "Test", x=analyte[0], y=[i]*len(analyte[0]), pen=None, symbolPen=pg.intColor(analyte_colour_array[colour],  values=5), symbolBrush=pg.intColor(analyte_colour_array[colour],values=5),symbol='o', symbolSize= 10)
								min_int = 50000
								int_x=40
								if max_int < min_int:
									test.setSymbolSize((max_int/230000)*int_x)
								else:
									test.setSymbolSize((max_int/230000)*4)
							else:
								pass
							test = self.ex_bOn_plot.plot(title= "Test", x=analyte[0], y=[i]*len(analyte[0]), pen=None, symbolPen=pg.intColor(analyte_colour_array[colour],  values=5, alpha = alpha), symbolBrush=pg.intColor(analyte_colour_array[colour],values=5, alpha = alpha),symbol='o', symbolSize= 10)
							min_int = 50000
							int_x=40
							if max_int < min_int:
								test.setSymbolSize((max_int/230000)*int_x)
							else:
								test.setSymbolSize((max_int/230000)*4)

							colour+=1

					


				
		
				# ticks = sample_names

				# print(ticks)
				# ticksdict = dict(enumerate(ticks))
				list1 = (sample_names)
				print(list1)


				print(list1)
				print(sample_names)

				ticks=[list(zip(range(len(sample_names)), list1))]

				# print('Diversity Test 3')
				ax = self.ex_bOn_plot.getAxis('left')

				ax.setTickSpacing(1,1)

				ax.setTicks(ticks)




			else:
				pass

			self.ex_bOn_plot.enableAutoRange(axis=None)
			print('PLOT EXPERIMENT END')
			self.c.append('Experiment Tab Loaded Succesfully')


			if len(sample_names) < 10:
				self.ex_bOn_plot.setYRange(0, -50)

			

	def Plot_Diversity(self):
		if open_state ==2:
			print('PLOT DIVERSITY BEGIN')
			print('Loading Diversity Tab')
			# self.div_plot.setYRange(min=-20,max=1)
			# self.div_plot.setXRange(min=0,max=140)
			self.div_plot.disableAutoRange(axis=None)
			tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
			Sample_State = tab_state.return_sample_state()
			Replicate_State = tab_state.return_replicate_state()
			Experiment_State = tab_state.return_experiment_state()
			Diversity_State = tab_state.return_diversity_state()


			if Diversity_State == False:
				comboAnalyte=self.comboAnalyte.currentText()

				print('Loading Diversity Tab')
				print('Importing DataFrame')
				sample_div, experiment_df = self.Diversity_Plot_DF()
				print('diversity samples', sample_div)
				self.div_plot.clear()
				
				for i in range(0, len(sample_div)):
					print('Loading Diversity Tab')
					analyte_array = []
					experiment_sample_df = experiment_df[experiment_df.sample_name == sample_div[i]]


					experiment_array_BOn = []

					print('Loading Diversity Tab')
					for experiment_analyte_id in experiment_sample_df.experiment_analyte_id.unique():
						experiment_array_BOn.append([experiment_sample_df[experiment_sample_df["experiment_analyte_id"] == experiment_analyte_id]["experiment_analyte_id"].to_numpy(), 
										experiment_sample_df[experiment_sample_df["experiment_analyte_id"] == experiment_analyte_id]["sample_analyte_max_intensity"].to_numpy()])
						analyte_array.append(experiment_analyte_id)
					
					
					analyte_array = [0 if str(z) == 'nan' else z for z in analyte_array]



					analyte_colour_array = [0]*len(analyte_array)

					print('Loading Diversity Tab')
					
					for a in range(0, len(analyte_array)):

						analyte_colour_array[a] = int((analyte_array[a] + 1) % 299)
						
					print('Dataframe Imported')



					if comboAnalyte == "Show All":


						print('Starting Plot')
						print('Loading Diversity Tab')
						colour = 0
						for analyte in experiment_array_BOn:

							print('Plotting Sample', [colour])




							y=[i]*len(analyte[0])

							x=analyte[0]

							max_int = max(analyte[1])


							test = self.div_plot.plot(title= "Test", x=x, y=y, pen=None, symbolPen=pg.intColor(analyte_colour_array[colour],  values=5), symbolBrush=pg.intColor(analyte_colour_array[colour],values=5),symbol='o')
							min_int = 50000
							int_x=40
							if max_int < min_int:
								test.setSymbolSize((max_int/230000)*int_x)
							else:
								test.setSymbolSize((max_int/230000)*4)

							colour+=1
						print('Plotting Finished')



					if comboAnalyte != "Show All":


						print('Loading Diversity Tab')
						alpha=30
						colour = 0

						for analyte in experiment_array_BOn:




							max_int = max(analyte[1])
							if analyte_array[colour] == float(comboAnalyte):


								test = self.div_plot.plot(title= "Test", x=analyte[0], y=[i]*len(analyte[0]), pen=None, symbolPen=pg.intColor(analyte_colour_array[colour],  values=5), symbolBrush=pg.intColor(analyte_colour_array[colour],values=5),symbol='o', symbolSize= 10)
								min_int = 50000
								int_x=40
								if max_int < min_int:
									test.setSymbolSize((max_int/230000)*int_x)
								else:
									test.setSymbolSize((max_int/230000)*4)
							else:
								pass
							test = self.div_plot.plot(title= "Test", x=analyte[0], y=[i]*len(analyte[0]), pen=None, symbolPen=pg.intColor(analyte_colour_array[colour],  values=5, alpha = alpha), symbolBrush=pg.intColor(analyte_colour_array[colour],values=5, alpha = alpha),symbol='o', symbolSize= 10)
							min_int = 50000
							int_x=40
							if max_int < min_int:
								test.setSymbolSize((max_int/230000)*int_x)
							else:
								test.setSymbolSize((max_int/230000)*4)

							colour+=1

					


				


				
				print('Diversity Test 1')
				# ticks = sample_names
				print('Diversity Test 2')
				# print(ticks)
				# ticksdict = dict(enumerate(ticks))
				list1 = (sample_names)
				print(list1)


				print(list1)
				print(sample_names)

				ticks=[list(zip(range(len(sample_names)), list1))]

				# print('Diversity Test 3')
				ax = self.div_plot.getAxis('left')
				print('Diversity Test 4')
				ax.setTickSpacing(1,1)
				print('Diversity Test 5')
				ax.setTicks(ticks)

				print('Diversity Test 6')
			else:
				pass
			print('Diversity Test 7')
			self.div_plot.enableAutoRange(axis=None)
			self.c.append('Diversity Tab Loaded Succesfully')
	##################################################################################################################	
			print('PLOT DIVERSITY END')
			if len(sample_names) < 10:
				self.div_plot.setYRange(0, -50)

		return 




	def Plot_ms1(self):
		if open_state ==2:
			print('PLOT MS1 BEGIN')
			chooseAnalyte = self.chooseAnalyte.currentText()
			print(chooseAnalyte)
			comboAnalyte=self.comboAnalyte.currentText()
			comboText=self.comboBox.currentText()

			print('PLOT MS1 BEGIN*')

			self.ms1_rightButton.setEnabled(False)
			self.ms1_leftButton.setEnabled(False)






			self.ms1_plot.clear()
			self.ms2_plot.clear()
			# self.ms1_plot.setLimits(ymin=0, ymax=100000)
			print('PLOT MS1 BEGIN**')
			df, df2 , df3= self.ms1_Plot_DF()

			# self.ms1_plot.disableAutoRange(axis=None)
			print('PLOT MS1 BEGIN**')
			if chooseAnalyte == "Replicate Analyte ID":
				self.ms1_plot.setXRange(min=0,max=1200)
				self.ms1_plot.setYRange(min=6, max=100)
				# self.ms1_plot.setLimits(ymin=0, ymax=100)
				ms1_data = []
				analyte_id_array = []

				if comboAnalyte != "Show All":
					
					for analyte_id in df.replicate_analyte_id.unique():
						ms1_data.append([df[df["replicate_analyte_id"] == analyte_id]["relative_intensity"].to_numpy(),
					                     df[df["replicate_analyte_id"] == analyte_id]["average_mass"].to_numpy()])
						analyte_id_array.append(analyte_id)
					analyte_colour_array = [0]*len(analyte_id_array)


					i = 0
					while i < len(analyte_id_array):
						analyte_colour_array[i] = int((analyte_id_array[i] + 1) % 299)
						i+=1



					colour=0
					self.ms1_plot.clear()
					for analyte in ms1_data:

						

						if float(comboAnalyte) == analyte_id_array[colour]:

							pen = pg.mkPen(color=(0, 0, 0))
							self.test = self.ms1_plot.plot( x=analyte[1], y=analyte[0], pen=pen )
							self.test.sigPointsClicked.connect(self.ms1_clicked)


							i = 0
							while i < len(analyte[1]):
								np.round(analyte[1],2)
								np.round(analyte[0],2)
								if analyte[0][i] > 10:





		
									text = pg.TextItem(str(np.round(analyte[1][i],4)),color=(0,0,0))
									self.ms1_plot.addItem(text)


									text.setPos(float(analyte[1][i]), float(analyte[0][i]) + 5)
									

									i+=1
								else:
									i+=1
						else:

							pass
						colour +=1

				# self.ms1_plot.enableAutoRange(axis=None)



			################# Sample Analyte ID ms1 Data #####################
			print('PLOT MS1 BEGIN***')
			if chooseAnalyte == "Sample Analyte ID":


				print('PLOT MS1 BEGIN****')

				ms1_data = []
				analyte_id_array = []

				if comboAnalyte != "Show All":





					df=df2[df2.analyte_id == float(comboAnalyte)]


					########## Set Max and Min #################

					mz_range_upper = df[df.mz == max(df.mz)].mz
					mz_range_upper = mz_range_upper.to_numpy()
					mz_range_lower = df[df.mz == min(df.mz)].mz
					mz_range_lower = mz_range_lower.to_numpy()

					self.ms1_plot.setXRange(min=mz_range_lower[0] - 5,max=mz_range_upper[0]+5)
					self.ms1_plot.setYRange(min=16000, max=300000)

					print('df sorted by chosen analyte',df)




					#### Find the max number in intensity column


					max_scan = df[df.intensity == df.intensity.max()].scan
					print('MAX SCAN DF',max_scan)


					max_scan_value = max_scan.to_numpy()
					print('MAX SCAN',max_scan_value)


					df= df[df.scan == max_scan_value[0]]
					print('DF sorted by max scan value',df)


					df.sort_values('scan', inplace=True)



					mz_array = df['mz'].to_numpy()

					print('mz array',mz_array)
					intensity_array = df['intensity'].to_numpy()
					print('intensity array',intensity_array)



					################# Create lists of zero's ###############

					intensity_zeros = [0]*len(intensity_array)
					mass_lower = [0]*len(mz_array)


					i=0
					while i < len(mz_array):
					    mass_lower[i] = mz_array[i] - 0.0001
					    i+=1

					mass_upper = [0]*len(mz_array)
					i=0
					while i < len(mz_array):
					    mass_upper[i] = mz_array[i] + 0.0001
					    i+=1




				    ################## Create Data Frames ####################




					##################### mz values ######################
					data1 = {'mz': mz_array,
					        'intensity': intensity_array
					        }

					##################### upper zeros ######################
					data2 = {'mz': mass_upper,
					        'intensity': intensity_zeros
					        }
					##################### lower zeros ######################
					data3 = {'mz': mass_lower,
					        'intensity': intensity_zeros
					        }

					df_1 = pandas.DataFrame (data1, columns = [ 'mz','intensity'])

					df_2 = pandas.DataFrame (data2, columns = ['mz','intensity'])

					df_3 = pandas.DataFrame (data3, columns = ['mz','intensity'])

					#### Combine Data Frames #######
					df_combine = [df_1,df_2,df_3]

					df_combine = pandas.concat(df_combine)

					df_combine.sort_values('mz', inplace=True)


					print("DF COMBINE 1", df_combine)





					mz_array = df_combine['mz'].to_numpy()
					intensity_array = df_combine['intensity'].to_numpy()
					for analyte_id in df2.analyte_id.unique():
						ms1_data.append([df2[df2["analyte_id"] == analyte_id]["intensity"].to_numpy(),
					                     df2[df2["analyte_id"] == analyte_id]["mz"].to_numpy()])
						analyte_id_array.append(analyte_id)

					

					analyte_colour_array = [0]*len(analyte_id_array)

					analyte_id_array = [0 if str(z) == 'nan' else z for z in analyte_id_array]
					i = 0
					while i < len(analyte_id_array):
						analyte_colour_array[i] = int((analyte_id_array[i] + 1) % 299)
						i+=1



					colour=0
					self.ms1_plot.clear()
					for analyte in ms1_data:
						np.round(mz_array[0],2)
						np.round(intensity_array[0],2)
						if float(comboAnalyte) == analyte_id_array[colour]:
				
							analyte[1] = [0 if str(i) == 'nan' else i for i in analyte[1]]
							analyte[0] = [0 if str(i) == 'nan' else i for i in analyte[0]]
							pen = pg.mkPen(color=(0, 0, 0))
							
							self.test = self.ms1_plot.plot( x=mz_array, y=intensity_array,pen=pen)




							i = 0
							while i < len(mz_array):
								if intensity_array[i] > 0:

									text1 = pg.TextItem(str(np.round(mz_array[i],4)),color=(0,0,0))
									self.ms1_plot.addItem(text1)


									text1.setPos(float(mz_array[i]), float(intensity_array[i])+15000)

									i+=1
								else:
									i+=1
						else:
							pass
						colour +=1


					with open((os.path.join(output_path, experiment_name + "_experiment_import_parameters.pickle")), 'rb') as f:
						experiment_info = pickle.load(f)

					if experiment_info.ms2_type == 'DDA':
						ms2 = import_ms2_dataframe(input_structure, input_type, comboText)

						ms2_sorted = ms2[ms2.analyte_id == float(comboAnalyte)]
						
						unique_mass = []

						for average_mass in ms2_sorted.ms1_average_mass.unique():
							unique_mass.append(average_mass)
						self.length_unique_mass = len(unique_mass)
						unique_mass = np.round_(unique_mass, decimals = 5)
						print('UNIQUE MASS',unique_mass)
						df_combine = df_combine.round(4)
						print('DF COMBINE',df_combine)
						test1, test2 , test3= self.ms1_Plot_DF()


						print('ms1 Data',test2[test2.analyte_id == 5])
						print('ms2',ms2[ms2.analyte_id == 5])




						unique_intensity = []
						j = 0





						if self.right_button_is_clicked > len(unique_mass) - 1:
							self.right_button_is_clicked = 0

						if len(unique_mass) == 1:
							self.right_button_is_clicked = 0



						print('BUTTON IS CLICKED',self.right_button_is_clicked)
						while j < len(unique_mass):




							mass_upper_lower = []
							
							intensity_upper_lower = []

							value = format(float(unique_mass[j]), '.4f')
							print(value)
							sample_analyte_df_sorted = df_combine[df_combine.mz == float(value)]



							mass_upper_lower.append(float(unique_mass[j]) - 0.0001)
							mass_upper_lower.append(float(unique_mass[j]))
							mass_upper_lower.append(float(unique_mass[j]) + 0.0001)

							print('SAMPLE ANLAYTE DF SORTED',sample_analyte_df_sorted)
							unique_intensity.append(sample_analyte_df_sorted.intensity.to_numpy())
							print('TEST1')
							intensity_upper_lower.append(0)
							print('TEST2')
							intensity_upper_lower.append(float(unique_intensity[j][0]))
							print('TEST3')
							intensity_upper_lower.append(0)
							print('TEST4')

							print(mass_upper_lower)
							print('TEST5')
							print(intensity_upper_lower)
							print('TEST6')
							# print("Unique Intensity", unique_intensity)
							# print("Unique Intensity", unique_intensity[0][0])
							# print(mass_upper_lower)
							# print(intensity_upper_lower)

							pen = pg.mkPen(color=(255, 0, 0),width = 1.5)
							print('TEST7')
							self.ms1_plot.plot( x=mass_upper_lower, y=intensity_upper_lower,pen=pen)
							print('TEST8')
							# button_is_clicked = repr(self.sender())

			

							j+=1

						mass_upper_lower = []
						print('TEST9')
						intensity_upper_lower = []
						print('TEST10')
						value = format(float(unique_mass[self.right_button_is_clicked]), '.4f')
						print('TEST11')
						print(value)
						print('TEST12')
						sample_analyte_df_sorted = df_combine[df_combine.mz == float(value)]
						print('TEST13')

						mass_upper_lower.append(float(unique_mass[self.right_button_is_clicked]) - 0.0001)
						print('TEST14')
						mass_upper_lower.append(float(unique_mass[self.right_button_is_clicked]))
						print('TEST15')

						mass_upper_lower.append(float(unique_mass[self.right_button_is_clicked]) + 0.0001)
						print('TEST16')



						unique_intensity.append(sample_analyte_df_sorted.intensity.to_numpy())
						print('TEST17')

						intensity_upper_lower.append(0)
						print('TEST18')
						# print('UNIQUE Intensity',unique_intensity[2][0])
						intensity_upper_lower.append(float(unique_intensity[self.right_button_is_clicked][0]))
						print('TEST19')

						intensity_upper_lower.append(0)
						print('TEST20')

						pen = pg.mkPen(color=(152,251,152),width = 1.5)
						print('TEST21')

						
						self.ms1_plot.plot( x=mass_upper_lower, y=intensity_upper_lower,pen=pen)


						if self.right_button_is_clicked == len(unique_mass) - 1:
							self.ms1_rightButton.setEnabled(False)

						else:
							self.ms1_rightButton.setEnabled(True)

						if self.right_button_is_clicked == 0:
							self.ms1_leftButton.setEnabled(False)

						else:
							self.ms1_leftButton.setEnabled(True)


						print('TEST22')
						print(mass_upper_lower)


						ms2_mass_data = []
						ms2_sorted=ms2_sorted[ms2_sorted.ms1_average_mass == float(unique_mass[self.right_button_is_clicked])]
						ms2_sorted.sort_values('ms2_data', inplace=True)
						max_ms2_mass = ms2_sorted.ms2_data.max()
						print(max_ms2_mass)
						min_ms2_mass = ms2_sorted.ms2_data.min()
						print(max_ms2_mass)
						print(ms2_sorted)

						ms2_mass_data = ms2_sorted.ms2_data.to_numpy()
						print(ms2_mass_data)

						ms2_intensity_data = ms2_sorted.Intensity.to_numpy()
						print(ms2_intensity_data)
						pen = pg.mkPen(color=(0, 0, 0))
						self.ms2_plot.setXRange(min=float(min_ms2_mass) - 0.1,max=float(max_ms2_mass) + 0.1)
						self.ms2_plot.setYRange(min=1200, max=20000)
						self.ms2_plot.clear()
						self.ms2_plot.plot( x=ms2_mass_data, y=ms2_intensity_data,pen=pen)

						print(ms2_mass_data)

						i = 0
						while i < len(ms2_mass_data):


							if ms2_intensity_data[i] > 0:

								text2 = pg.TextItem(str(np.round(ms2_mass_data[i],4)),color=(0,0,0))
								self.ms2_plot.addItem(text2)


								text2.setPos(float(ms2_mass_data[i]), float(ms2_intensity_data[i])+1500)

								i+=1
							else:
								i+=1

			if chooseAnalyte == "Experiment Analyte ID":
				self.ms1_plot.setXRange(min=0,max=1200)
				self.ms1_plot.setYRange(min=5, max=100)
				ms1_data = []
				analyte_id_array = []

				if comboAnalyte != "Show All":

					for analyte_id in df3.experiment_analyte_id.unique():
						ms1_data.append([df3[df3["experiment_analyte_id"] == analyte_id]["relative_intensity"].to_numpy(),
					                     df3[df3["experiment_analyte_id"] == analyte_id]["average_mass"].to_numpy()])
						analyte_id_array.append(analyte_id)
					analyte_colour_array = [0]*len(analyte_id_array)


					i = 0
					while i < len(analyte_id_array):
						analyte_colour_array[i] = int((analyte_id_array[i] + 1) % 299)
						i+=1



					colour=0
					self.ms1_plot.clear()
					for analyte in ms1_data:

						if float(comboAnalyte) == analyte_id_array[colour]:

							pen = pg.mkPen(color=(0, 0, 0))
							self.test = self.ms1_plot.plot( x=analyte[1], y=analyte[0],pen=pen)
							self.test.sigPointsClicked.connect(self.ms1_clicked)
							i = 0
							while i < len(analyte[1]):
								if analyte[0][i] > 0:

									text = pg.TextItem(str(np.round(analyte[1][i],4)),color=(0,0,0))
									self.ms1_plot.addItem(text)


									text.setPos(float(analyte[1][i]), float(analyte[0][i])+0.1)
									i+=1
								else:
									i+=1
						else:
							pass

						colour +=1
					
				# self.ms1_plot.enableAutoRange(axis=None)
			print('PLOT MS1 END')
		return











	def Plot_ms2(self):

		if open_state ==2:

			with open((os.path.join(output_path, experiment_name + "_experiment_import_parameters.pickle")), 'rb') as f:
				experiment_info = pickle.load(f)




				print('PLOT MS2 BEGIN')
				chooseAnalyte = self.chooseAnalyte.currentText()
				print('FILL COMBO1')
				comboAnalyte=self.comboAnalyte.currentText()
				print('FILL COMBO2')
				comboText=self.comboBox.currentText()
				print('FILL COMBO3')

				if self.lock_aspect.isChecked():
					print('FILL COMBO4')
					self.ms2_plot.setXLink(self.ms1_plot)
					self.lock_aspect.setIcon(QIcon('Lock_icon.png'))
					print('FILL COMBO5')
				# if it is unchecked 

				else: 
					print('FILL COMBO6')
					self.ms2_plot.setXLink(self.ms2_plot)
					self.lock_aspect.setIcon(QIcon('Unlock_icon.png'))
					print('FILL COMBO7')

			print('FILL COMBO8')
			print(experiment_info.ms2_type)
			if experiment_info.ms2_type == 'DIA':
		




				print('FILL COMBO?')


				self.ms2_plot.clear()


				# self.ms1_plot.setLimits(ymin=0, ymax=100000)
				df= self.ms2_Plot_DF()
				print('DATA FRAME ms2',df)
				# self.ms1_plot.disableAutoRange(axis=None)
				print('PLOT MS2 BEGIN')
				if chooseAnalyte == "Replicate Analyte ID":
					print('PLOT MS2 BEGIN')
					self.ms2_plot.setXRange(min=0,max=1200)
					self.ms2_plot.setYRange(min=6, max=100)
					# self.ms1_plot.setLimits(ymin=0, ymax=100)
					ms2_data = []
					analyte_id_array = []
					print('PLOT MS2 BEGIN')
					if comboAnalyte != "Show All":
						print('PLOT MS2 BEGIN')

						for analyte_id in df.replicate_analyte_id.unique():

							ms2_data.append([df[df["replicate_analyte_id"] == analyte_id]["relative_intensity"].to_numpy(),
						                     df[df["replicate_analyte_id"] == analyte_id]["average_mass"].to_numpy()])
							analyte_id_array.append(analyte_id)
						analyte_colour_array = [0]*len(analyte_id_array)
						print('PLOT MS2 BEGIN')


						i = 0
						while i < len(analyte_id_array):
							analyte_colour_array[i] = int((analyte_id_array[i] + 1) % 299)
							i+=1
						print('PLOT MS2 BEGIN')


						colour=0
						self.ms2_plot.clear()
						print('PLOT MS2 BEGIN')
						for analyte in ms2_data:
							print('PLOT MS2 BEGIN')
							

							if float(comboAnalyte) == analyte_id_array[colour]:

								pen = pg.mkPen(color=(0, 0, 0))
								self.test = self.ms2_plot.plot( x=analyte[1], y=analyte[0], pen=pen )



								i = 0
								while i < len(analyte[1]):
									np.round(analyte[1],2)
									np.round(analyte[0],2)
									if analyte[0][i] > 0:
										text = pg.TextItem(str(np.round(analyte[1][i],4)),color=(0,0,0))
										self.ms2_plot.addItem(text)


										text.setPos(float(analyte[1][i]), float(analyte[0][i]) + 5)
										i+=1
									else:
										i+=1
							else:

								pass
							colour +=1

					# self.ms1_plot.enableAutoRange(axis=None)




				if chooseAnalyte == "Sample Analyte ID":
					pass




				if chooseAnalyte == "Experiment Analyte ID":
					pass
						
					# self.ms1_plot.enableAutoRange(axis=None)
				print('PLOT MS2 END')



			if experiment_info.ms2_type == 'DDA':
				print('DDA')
		return











	def Plot_Analyte_Legend(self):


		if open_state ==2:
			print('PLOT ANALYTE LEGEND BEGIN')
			self.analyte_legend.disableAutoRange(axis=None)
			tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
			Sample_State = tab_state.return_sample_state()
			Replicate_State = tab_state.return_replicate_state()
			Experiment_State = tab_state.return_experiment_state()
			Diversity_State = tab_state.return_diversity_state()

			chooseAnalyte = self.chooseAnalyte.currentText()
			comboText=self.comboBox.currentText()

			sorted_df_mz = import_dataframe(output_path, input_type, comboText)

			# sorted_df_mz, sorted_df_rt= self.Sample_Plot_DF()


			mz_array = []
			analyte_mz_array = []

			self.analyte_legend.clear()
			if chooseAnalyte == "Sample Analyte ID":
				print('Sample Analyte')

				for analyte_id in sorted_df_mz.analyte_id.unique():

					mz_array.append([sorted_df_mz[sorted_df_mz["analyte_id"] == analyte_id]["analyte_id"].to_numpy()])

					analyte_mz_array.append(analyte_id)



				analyte_mz_array = [0 if str(i) == 'nan' else i for i in analyte_mz_array]





				analyte_colour_mz_array = [0]*len(analyte_mz_array)


				i = 0
				while i < len(analyte_mz_array):
					analyte_colour_mz_array[i] = int((analyte_mz_array[i] + 1) % 299)
					i+=1

			if chooseAnalyte == "Replicate Analyte ID":
				print('Replicate Analyte')

				for analyte_id in sorted_df_mz.analyte_id.unique():

					mz_array.append([sorted_df_mz[sorted_df_mz["replicate_analyte_id"] == analyte_id]["replicate_analyte_id"].to_numpy()])

					analyte_mz_array.append(analyte_id)



				analyte_mz_array = [0 if str(i) == 'nan' else i for i in analyte_mz_array]





				analyte_colour_mz_array = [0]*len(analyte_mz_array)


				i = 0
				while i < len(analyte_mz_array):
					analyte_colour_mz_array[i] = int((analyte_mz_array[i] + 1) % 299)
					i+=1


			if chooseAnalyte == "Experiment Analyte ID":
				print('Experiment Analyte')

				for analyte_id in sorted_df_mz.analyte_id.unique():

					mz_array.append([sorted_df_mz[sorted_df_mz["experiment_analyte_id"] == analyte_id]["experiment_analyte_id"].to_numpy()])

					analyte_mz_array.append(analyte_id)



				analyte_mz_array = [0 if str(i) == 'nan' else i for i in analyte_mz_array]





				analyte_colour_mz_array = [0]*len(analyte_mz_array)


				i = 0
				while i < len(analyte_mz_array):
					analyte_colour_mz_array[i] = int((analyte_mz_array[i] + 1) % 299)
					i+=1



			colour=0

			for analyte in analyte_mz_array:


				#self.mz_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_array[colour],  values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_array[colour],values=5,alpha=alpha),symbol='o', symbolSize=3)
				#colour+=1

				x=[0]
				y=[float(analyte),0]
				


				self.test = self.analyte_legend.plot( x=x, y=y, pen=None, symbolPen=pg.intColor( analyte_colour_mz_array[colour],  values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_mz_array[colour],values=5,alpha=255),symbol='o', symbolSize=6)
				self.test.sigPointsClicked.connect(self.clicked_Legend)
				colour+=1


				text = pg.TextItem(str(int(y[0])),color=(0,0,0))

				self.analyte_legend.addItem(text)


				text.setPos(0.1, y[0]-0.70)

			self.analyte_legend.getViewBox().invertY(True)
			self.analyte_legend.setLimits(xMin=-.05, xMax=0.3,yMin=0.5, yMax=max(analyte_mz_array)+20)

			print('PLOT ANALYTE LEGEND END')
			self.analyte_legend.enableAutoRange(axis=None)

			self.analyte_legend.setYRange(min=2, max=40)
			self.analyte_legend.setXRange(min=0, max=0.5)
			self.analyte_legend.setMouseEnabled(x=False, y=True)

	def reset_all(self):
			state = State(self.Nullbutton.setChecked(False), self.Blankbutton.setChecked(False), 0, 1200, 0, 7,self.Toggle_rt.setChecked(False))
			self.Blankbutton.setStyleSheet("background-color : lightgreen") 
			self.Blankbutton.setText('On')
			self.Nullbutton.setStyleSheet("background-color : lightgreen") 
			self.Nullbutton.setText('On')
			self.Toggle_rt.setText('Toggle Scatter Plot')
			self.Toggle_rt.setStyleSheet("background-color : lightgrey")
			self.massCheckbox.setChecked(False)

			self.mMinBox.setText('0')
			self.mMaxBox.setText('1200')
			self.rtMinBox.setText('0')
			self.rtMaxBox.setText('7')



	def reset_massRange(self):

		state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), 0, 1200, self.rtMinBox.text(), self.rtMaxBox.text(),self.Toggle_rt.isChecked())
		self.mMinBox.setText('0')
		self.mMaxBox.setText('1200')

	def reset_rtRange(self):
		state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), 0, 7,self.Toggle_rt.isChecked())
		self.rtMinBox.setText('0')
		self.rtMaxBox.setText('7')


	def Sample_Plot_DF(self):





		if open_state == 2:

			print(sample_names)
			comboText=self.comboBox.currentText()
			print(self.comboBox)


			print(self.comboBox.currentText())
			print('SAMPLE PLOT DF BEGIN')
			chooseAnalyte = self.chooseAnalyte.currentText()
			tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
			Sample_State = tab_state.return_sample_state()
			Replicate_State = tab_state.return_replicate_state()
			Experiment_State = tab_state.return_experiment_state()



			print('TEST1')

			print(input_structure)
			print('TEST2')
			print(input_type)
			print('TEST3')
			print(comboText)
			print('TEST4')
			print(output_path)
			total_df_test = import_dataframe(output_path, input_type, comboText)
			print('TEST2')

			max_mass = total_df_test[total_df_test.mz == total_df_test.mz.max()].mz
			global max_mass_value
			max_mass_value = max_mass.to_numpy()

			max_rt = total_df_test[total_df_test.rt == total_df_test.rt.max()].rt
			global max_rt_value
			max_rt_value = max_rt.to_numpy()
			if Sample_State == False:
				if self.massCheckbox.isChecked() == False:
					# if self.rtCheckbox.isChecked() == False:
					state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), 0, max_mass_value[0], 0, max_rt_value[0],self.Toggle_rt.isChecked())

				#if self.rtCheckbox.isChecked() == False:
				#	state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), 0, 7)

				if self.massCheckbox.isChecked() == True:
					# if self.rtCheckbox.isChecked() == False:
					state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), self.rtMinBox.text(), self.rtMaxBox.text(),self.Toggle_rt.isChecked())



				# if self.rtCheckbox.isChecked() == True:
				# 	if self.massCheckbox.isChecked() == False:
				# 		state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), self.rtMinBox.text(), self.rtMaxBox.text())




				BlankState = state.return_blank_state()
				NullState = state.return_null_state()
				Toggle_rt = state.return_Toggle_rt()
				comboText=self.comboBox.currentText()

				############### Fetch data for mz_plot and plot sample from comboBox ###############

				total_df = import_dataframe(output_path, input_type, comboText)
				print('ANALYTE ID 3',total_df[total_df.analyte_id == 2])

				print(chooseAnalyte)

				if chooseAnalyte == "Sample Analyte ID":
					print('Blank Test 1')
					if NullState == True:

						total_df = total_df.dropna(subset=['analyte_id'])
						total_df=total_df[total_df.analyte_id != None]

						
					else:
						pass


					if BlankState == True:
						print('Blank Test 2')
						total_df=total_df[total_df.blank_analyte != False]
						print('Blank Test 3')
					else:
						pass



					print('Blank Test 4')


					total_df=total_df[total_df.replicate != 2]
					total_df=total_df[total_df.replicate != 3]
					print('Blank Test 3')
					total_mz_df = total_df[total_df.rt > float(state.return_rtMin_state())]
					total_mz_df = total_mz_df[total_mz_df.rt < float(state.return_rtMax_state())]

					sorted_df_mz = total_mz_df.sort_values(by=['mz'])
					print('SORTED DF MZ',sorted_df_mz[sorted_df_mz.analyte_id == 2])

					total_rt_df = total_df[total_df.mz > float(state.return_massMin_state())]

					total_rt_df = total_rt_df[total_rt_df.mz < float(state.return_massMax_state())]
					sorted_df_rt = total_rt_df.sort_values(by=['rt'])
			




				if chooseAnalyte == "Replicate Analyte ID":
					if NullState == True:
						total_df = total_df.dropna(subset=['replicate_analyte_id'])
						total_df=total_df[total_df.replicate_analyte_id != None]

					else:
						pass
					if BlankState == True:
						total_df=total_df[total_df.blank_analyte != True]
					else:
						pass
					total_df=total_df[total_df.replicate != 2]
					total_df=total_df[total_df.replicate != 3]

					total_mz_df = total_df[total_df.rt > float(state.return_rtMin_state())]
					total_mz_df = total_mz_df[total_mz_df.rt < float(state.return_rtMax_state())]

					sorted_df_mz = total_mz_df.sort_values(by=['mz'])


					total_rt_df = total_df[total_df.mz > float(state.return_massMin_state())]

					total_rt_df = total_rt_df[total_rt_df.mz < float(state.return_massMax_state())]
					sorted_df_rt = total_rt_df.sort_values(by=['rt'])



				if chooseAnalyte == "Experiment Analyte ID":
					if NullState == True:
						total_df = total_df.dropna(subset=['experiment_analyte_id'])
						total_df=total_df[total_df.experiment_analyte_id != None]

					else:
						pass
					if BlankState == True:
						total_df=total_df[total_df.blank_analyte != True]
					else:
						pass
					total_df=total_df[total_df.replicate != 2]
					total_df=total_df[total_df.replicate != 3]

					total_mz_df = total_df[total_df.rt > float(state.return_rtMin_state())]
					total_mz_df = total_mz_df[total_mz_df.rt < float(state.return_rtMax_state())]

					sorted_df_mz = total_mz_df.sort_values(by=['mz'])


					total_rt_df = total_df[total_df.mz > float(state.return_massMin_state())]

					total_rt_df = total_rt_df[total_rt_df.mz < float(state.return_massMax_state())]
					sorted_df_rt = total_rt_df.sort_values(by=['rt'])
				# else:
				# 	pass
			else:
				sorted_df_mz = []
				sorted_df_rt = []




		

				# max_scan_value = max_scan.to_numpy()

				# df=df[df.scan == max_scan_value[0]]
			print('SAMPLE PLOT DF END')


		else:
			sorted_df_mz = []

			sorted_df_rt = []


		return sorted_df_mz, sorted_df_rt



	def Replicate_Plot_DF(self):


		if open_state == 2:


			print('REPLICATE PLOT DF BEGIN')
			chooseAnalyte = self.chooseAnalyte.currentText()
			tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
			Sample_State = tab_state.return_sample_state()
			Replicate_State = tab_state.return_replicate_state()
			Experiment_State = tab_state.return_experiment_state()


			comboText=self.comboBox.currentText()



			if Replicate_State == False:
				if self.massCheckbox.isChecked() == False:
					# if self.rtCheckbox.isChecked() == False:
					state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), 0, float(max_mass_value[0]), 0,float(max_rt_value[0]),self.Toggle_rt.isChecked())

				#if self.rtCheckbox.isChecked() == False:
				#	state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), 0, 7)

				if self.massCheckbox.isChecked() == True:
					# if self.rtCheckbox.isChecked() == False:
					state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), self.rtMinBox.text(), self.rtMaxBox.text(),self.Toggle_rt.isChecked())



				# if self.rtCheckbox.isChecked() == True:
				# 	if self.massCheckbox.isChecked() == False:
				# 		state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), self.rtMinBox.text(), self.rtMaxBox.text())




				BlankState = state.return_blank_state()
				NullState = state.return_null_state()
				comboText=self.comboBox.currentText()

				############### Fetch data for mz_plot and plot sample from comboBox ###############

				total_df = import_dataframe(output_path, input_type, comboText)

				if chooseAnalyte == "Sample Analyte ID":
					if NullState == True:
						total_df = total_df.dropna(subset=['analyte_id'])
						total_df=total_df[total_df.analyte_id != None]
					else:
						pass
					if BlankState == True:
						total_df=total_df[total_df.blank_analyte != True]
					else:
						pass
					total_df_r1=total_df[total_df.replicate != 2]
					total_df_r1=total_df_r1[total_df_r1.replicate != 3]

					total_mz_df_r1 = total_df_r1[total_df_r1.rt > float(state.return_rtMin_state())]
					total_mz_df_r1 = total_mz_df_r1[total_mz_df_r1.rt < float(state.return_rtMax_state())]



					sorted_df_mz_r1 = total_mz_df_r1.sort_values(by=['mz'])


					total_rt_df_r1 = total_df_r1[total_df_r1.mz > float(state.return_massMin_state())]

					total_rt_df_r1 = total_rt_df_r1[total_rt_df_r1.mz < float(state.return_massMax_state())]
					sorted_df_rt_r1 = total_rt_df_r1.sort_values(by=['rt'])

					total_df_r2=total_df[total_df.replicate != 1]
					total_df_r2=total_df_r2[total_df_r2.replicate != 3]

					total_mz_df_r2 = total_df_r2[total_df_r2.rt > float(state.return_rtMin_state())]
					total_mz_df_r2 = total_mz_df_r2[total_mz_df_r2.rt < float(state.return_rtMax_state())]



					sorted_df_mz_r2 = total_mz_df_r2.sort_values(by=['mz'])


					total_rt_df_r2 = total_df_r2[total_df_r2.mz > float(state.return_massMin_state())]

					total_rt_df_r2 = total_rt_df_r2[total_rt_df_r2.mz < float(state.return_massMax_state())]
					sorted_df_rt_r2 = total_rt_df_r2.sort_values(by=['rt'])




					total_df_r3=total_df[total_df.replicate != 1]
					total_df_r3=total_df_r3[total_df_r3.replicate != 2]

					total_mz_df_r3 = total_df_r3[total_df_r3.rt > float(state.return_rtMin_state())]
					total_mz_df_r3 = total_mz_df_r3[total_mz_df_r3.rt < float(state.return_rtMax_state())]



					sorted_df_mz_r3 = total_mz_df_r3.sort_values(by=['mz'])


					total_rt_df_r3 = total_df_r3[total_df_r3.mz > float(state.return_massMin_state())]

					total_rt_df_r3 = total_rt_df_r3[total_rt_df_r3.mz < float(state.return_massMax_state())]
					sorted_df_rt_r3 = total_rt_df_r3.sort_values(by=['rt'])
					#sorted_df_rt = total_df.sort_values(by=['rt'])

				if chooseAnalyte == "Replicate Analyte ID":
					if NullState == True:
						total_df = total_df.dropna(subset=['replicate_analyte_id'])
						total_df=total_df[total_df.replicate_analyte_id != None]
					else:
						pass
					if BlankState == True:
						total_df=total_df[total_df.blank_analyte != True]
					else:
						pass
					total_df_r1=total_df[total_df.replicate != 2]
					total_df_r1=total_df_r1[total_df_r1.replicate != 3]

					total_mz_df_r1 = total_df_r1[total_df_r1.rt > float(state.return_rtMin_state())]
					total_mz_df_r1 = total_mz_df_r1[total_mz_df_r1.rt < float(state.return_rtMax_state())]



					sorted_df_mz_r1 = total_mz_df_r1.sort_values(by=['mz'])


					total_rt_df_r1 = total_df_r1[total_df_r1.mz > float(state.return_massMin_state())]

					total_rt_df_r1 = total_rt_df_r1[total_rt_df_r1.mz < float(state.return_massMax_state())]
					sorted_df_rt_r1 = total_rt_df_r1.sort_values(by=['rt'])

					total_df_r2=total_df[total_df.replicate != 1]
					total_df_r2=total_df_r2[total_df_r2.replicate != 3]

					total_mz_df_r2 = total_df_r2[total_df_r2.rt > float(state.return_rtMin_state())]
					total_mz_df_r2 = total_mz_df_r2[total_mz_df_r2.rt < float(state.return_rtMax_state())]



					sorted_df_mz_r2 = total_mz_df_r2.sort_values(by=['mz'])


					total_rt_df_r2 = total_df_r2[total_df_r2.mz > float(state.return_massMin_state())]

					total_rt_df_r2 = total_rt_df_r2[total_rt_df_r2.mz < float(state.return_massMax_state())]
					sorted_df_rt_r2 = total_rt_df_r2.sort_values(by=['rt'])




					total_df_r3=total_df[total_df.replicate != 1]
					total_df_r3=total_df_r3[total_df_r3.replicate != 2]

					total_mz_df_r3 = total_df_r3[total_df_r3.rt > float(state.return_rtMin_state())]
					total_mz_df_r3 = total_mz_df_r3[total_mz_df_r3.rt < float(state.return_rtMax_state())]



					sorted_df_mz_r3 = total_mz_df_r3.sort_values(by=['mz'])


					total_rt_df_r3 = total_df_r3[total_df_r3.mz > float(state.return_massMin_state())]

					total_rt_df_r3 = total_rt_df_r3[total_rt_df_r3.mz < float(state.return_massMax_state())]
					sorted_df_rt_r3 = total_rt_df_r3.sort_values(by=['rt'])

					#sorted_df_rt = total_df.sort_values(by=['rt'])


				if chooseAnalyte== "Experiment Analyte ID":
					if NullState == True:
						total_df = total_df.dropna(subset=['experiment_analyte_id'])
						total_df=total_df[total_df.experiment_analyte_id != None]
					else:
						pass
					if BlankState == True:
						total_df=total_df[total_df.blank_analyte != True]
					else:
						pass
					total_df_r1=total_df[total_df.replicate != 2]
					total_df_r1=total_df_r1[total_df_r1.replicate != 3]

					total_mz_df_r1 = total_df_r1[total_df_r1.rt > float(state.return_rtMin_state())]
					total_mz_df_r1 = total_mz_df_r1[total_mz_df_r1.rt < float(state.return_rtMax_state())]



					sorted_df_mz_r1 = total_mz_df_r1.sort_values(by=['mz'])


					total_rt_df_r1 = total_df_r1[total_df_r1.mz > float(state.return_massMin_state())]

					total_rt_df_r1 = total_rt_df_r1[total_rt_df_r1.mz < float(state.return_massMax_state())]
					sorted_df_rt_r1 = total_rt_df_r1.sort_values(by=['rt'])

					total_df_r2=total_df[total_df.replicate != 1]
					total_df_r2=total_df_r2[total_df_r2.replicate != 3]

					total_mz_df_r2 = total_df_r2[total_df_r2.rt > float(state.return_rtMin_state())]
					total_mz_df_r2 = total_mz_df_r2[total_mz_df_r2.rt < float(state.return_rtMax_state())]



					sorted_df_mz_r2 = total_mz_df_r2.sort_values(by=['mz'])


					total_rt_df_r2 = total_df_r2[total_df_r2.mz > float(state.return_massMin_state())]

					total_rt_df_r2 = total_rt_df_r2[total_rt_df_r2.mz < float(state.return_massMax_state())]
					sorted_df_rt_r2 = total_rt_df_r2.sort_values(by=['rt'])




					total_df_r3=total_df[total_df.replicate != 1]
					total_df_r3=total_df_r3[total_df_r3.replicate != 2]

					total_mz_df_r3 = total_df_r3[total_df_r3.rt > float(state.return_rtMin_state())]
					total_mz_df_r3 = total_mz_df_r3[total_mz_df_r3.rt < float(state.return_rtMax_state())]



					sorted_df_mz_r3 = total_mz_df_r3.sort_values(by=['mz'])


					total_rt_df_r3 = total_df_r3[total_df_r3.mz > float(state.return_massMin_state())]

					total_rt_df_r3 = total_rt_df_r3[total_rt_df_r3.mz < float(state.return_massMax_state())]
					sorted_df_rt_r3 = total_rt_df_r3.sort_values(by=['rt'])
			else:
				sorted_df_mz_r1 = []
				sorted_df_mz_r2 = []
				sorted_df_mz_r3 = []
				sorted_df_rt_r1 = []
				sorted_df_rt_r2 = []
				sorted_df_rt_r3 = []
			print('REPLICATE PLOT DF END')
		else:
			sorted_df_mz_r1 = []
			sorted_df_rt_r1 = []
			sorted_df_mz_r2 = []
			sorted_df_rt_r2 = []
			sorted_df_mz_r3 = []
			sorted_df_rt_r3 = []
		print(sorted_df_mz_r1)
		# test_df = sorted_df_mz_r1[sorted_df_mz_r1.analyte_id == 1]
		# print(test_df)
		return sorted_df_mz_r1, sorted_df_rt_r1, sorted_df_mz_r2, sorted_df_rt_r2, sorted_df_mz_r3, sorted_df_rt_r3
################### Connecting Region 1 to Region 3 ###########################
	def Experiment_Plot_DF(self):


		if open_state ==2:


			print('EXPERIMENT PLOT DF BEGIN')
			tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
			Sample_State = tab_state.return_sample_state()
			Replicate_State = tab_state.return_replicate_state()
			Experiment_State = tab_state.return_experiment_state()
			

			if Experiment_State == False:
				print('EXPERIMENT PLOT DF 2')
				state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), self.rtMinBox.text(), self.rtMaxBox.text(),self.Toggle_rt.isChecked())
				print('Experiment Test 1')
				BlankState = state.return_blank_state()
				print('Experiment Test 2')
				NullState = state.return_null_state()
				print('Experiment Test 3')

				experiment_df = import_experiment_dataframe(output_path)
				print('Experiment Test 4')
				experiment_df =experiment_df.sort_values(by='sample_processing_order')
				print('Experiment Test 5')
				
				if BlankState == True:
					experiment_df=experiment_df[experiment_df.experiment_analyte_is_blank != True]
				else:
					pass
				sample_name_df = experiment_df.sample_name
				sample_name_sorted_df = sample_name_df.drop_duplicates()

				samples = sample_name_sorted_df.values
			else:
				samples=[]
				experiment_df = []
			print('EXPERIMENT PLOT DF END')



		else:

			samples = []

			experiment_df = []
		print(experiment_df)

		return samples, experiment_df




	def Diversity_Plot_DF(self):



		if open_state ==2:


			tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
			Sample_State = tab_state.return_sample_state()
			Replicate_State = tab_state.return_replicate_state()
			Experiment_State = tab_state.return_experiment_state()
			Diversity_State = tab_state.return_diversity_state()

			if Diversity_State == False:

				state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), self.rtMinBox.text(), self.rtMaxBox.text(),self.Toggle_rt.isChecked())

				BlankState = state.return_blank_state()
				NullState = state.return_null_state()

				experiment_df = import_experiment_dataframe(output_path)
				experiment_df =experiment_df.sort_values(by='sample_processing_order')
				
				if BlankState == True:
					experiment_df=experiment_df[experiment_df.experiment_analyte_is_blank != True]
				else:
					pass
				sample_name_df = experiment_df.sample_name
				sample_name_sorted_df = sample_name_df.drop_duplicates()
				sample_name_sorted_df = sample_name_sorted_df.iloc[::-1]

				samples = sample_name_sorted_df.values

			else:
				samples=[]
				experiment_df = []

			return samples, experiment_df

	def ms1_Plot_DF(self):


		if open_state == 2:
			print('MS1 PLOT DF')
			comboText=self.comboBox.currentText()
			print(comboText)
			print(output_path)
			df1, df2, df3 = import_ms1_dataframe(output_path, input_type, comboText)
			print('MS1 PLOT DF*')




			df4 = import_dataframe(output_path, input_type, comboText)   ############# Data frame for Sample level ms1 data

			print('MS1 PLOT DF**')



			df1.sort_values('average_mass', inplace=True)
			df3.sort_values('average_mass', inplace=True)



			################### include only first replicate for sample level ms1 data
			print('original dataframe analyte 5',df4[df4.analyte_id == 8])
			df=df4[df4.replicate == 1]




		else:

			df1 = []
			df = []
			df3 = []

		return df1, df, df3




	def is_right_button_clicked(self):
		print('***************************************************************UNIQUE MASS',self.length_unique_mass)

		if self.right_button_is_clicked != (self.length_unique_mass - 1):

			self.right_button_is_clicked +=1


		return self.right_button_is_clicked

	def is_left_button_clicked(self):

		if self.right_button_is_clicked != 0:
			self.right_button_is_clicked -=1

		else:
			self.right_button_is_clicked =0

		return self.right_button_is_clicked








	def openFileNamesDialog(self):
		global input_type
		input_type = "Samples"
		print(input_type)
		global open_state
		flags = QtWidgets.QFileDialog.ShowDirsOnly
		self.dir_ = QtWidgets.QFileDialog.getExistingDirectory(
			None,
			"Select directory",
			"",
			flags)
		if self.dir_:

			
			open_state = 2


			# self.path_manager.set_output_dirname(self.dir_)
			# self.refresh_paths()



			self.path_name_list = []

			for path in Path(self.dir_).rglob('*parameters.pickle'):
				
				self.path_name_list.append(path.name)

			global output_path
			output_path = self.dir_
			self.output_path_name = str(self.path_name_list[0])
			print('dir_',self.dir_)
			print('output_path_name',self.output_path_name)



			self.object1 = pandas.read_pickle(os.path.join(self.dir_, self.output_path_name))
			print("TEST1")
			global input_structure
			input_structure = self.object1
			global experiment_name
			experiment_name = input_structure.experiment_name
			print(input_structure)
			print('input type 1', input_type)
			print('EXPERIMENT NAME',experiment_name)
			print('input type 2', input_type)
			self.dir_sample = self.dir_[:self.dir_.find('Output')].strip()
			self.dir_sample = self.dir_sample + 'Samples'
			print(self.dir_sample)


			sample_name_list = []

			if input_type == "Samples":
			    base_input_directory = self.dir_sample

			else:
				print("ERROR - data_import.name_extract: Replicate compare. Sample input type not valid. Must be 'Samples' or "
						"'Blanks'")
				sys.exit()

			# Create list of sample names to analyze
			for sample in glob.glob(os.path.join(base_input_directory, "*_R[0-9]." + input_structure.ms_data_file_suffix)):
				sample_name = str(os.path.basename(sample).rsplit("_", maxsplit=1)[0])
				if sample_name not in sample_name_list:
					sample_name_list.append(sample_name)




			print(sample_name_list)




			global sample_names
			sample_names = sample_name_list
			print(self.object1.sample_directory)
			# print(self.input_type)
			print('SAMPLE NAMES',sample_names)

			# object1 = pandas.read_pickle((os.path.join(data_import.input_data_structure().output_directory, path_name)))
			n = 0
			print('COMBOBOX TEST1')
			while n < len(sample_names):
				print('COMBOBOX TEST2')
				print(sample_names[n])
				print('CURRENT TEXT',self.comboBox.currentText())
				# self.comboBox.addItem(str(sample_names[n]))
				print('COMBOBOX TEST3')
				
				n+=1
				print('COMBOBOX TEST4')
			print('CURRENT TEXT',self.comboBox.currentText())
		else:

			open_state = 1
		# with open((os.path.join(output_path, experiment_name + "_experiment_import_parameters.pickle")), 'rb') as f:
		# 	experiment_info = pickle.load(f)
		# if experiment_info.ms2_type == 'DDA':
		# 	self.mgf_action.setDisabled(True)

		# 	self.export_mgf.setDisabled(True)

		# else:
		# 	self.mgf_action.setDisabled(False)
		# 	self.export_mgf.setDisabled(False)

		print(input_structure.blanks_exist)
		if input_structure.blanks_exist == False:

			self.Blankbutton.setDisabled(True)

		else:
			self.Blankbutton.setDisabled(False)

	def info_test(self):


		if open_state == 2:

			print('TEST comboBox 1')
			print(sample_names)
			print('TEST comboBox 2')
			print(self.comboBox.currentText())
			print('TEST comboBox 3')
			
			print('TEST comboBox 4')
			print(self.comboBox.currentText())
			print('TEST comboBox 5')

			self.comboBox.blockSignals(True)
			self.comboBox.clear()
			self.comboBox.addItems(sample_names)
			self.comboBox.blockSignals(False)

			# n=0

			# while n < len(sample_names):
			# 	self.comboBox.addItem(str(sample_names[n]))

			# 	n+=1


			comboText = self.comboBox.currentText()
			total_df = import_dataframe(output_path, input_type, comboText)






			sorted_df = total_df.dropna(subset=['analyte_id'])
			sorted_df=total_df[total_df.analyte_id != None]
			sorted_df = sorted_df.sort_values(by=['analyte_id'])

			mz_array=[]
			analyte_array=[]

			for analyte_id in sorted_df.analyte_id.unique():


				analyte_array.append(int(analyte_id))



			n = 0



			print('ANALYTE ARRAY',analyte_array)

			self.comboAnalyte.blockSignals(True)
			self.comboAnalyte.clear()

			self.comboAnalyte.addItem('Show All')
			while n < len(analyte_array):

				self.comboAnalyte.addItem(str(analyte_array[n]))
				
				n+=1
			self.comboAnalyte.blockSignals(False)

			print('test test test')




			print('===== Test Total DF')
			total_df_test = import_dataframe(output_path, input_type, comboText)


			max_mass = total_df_test[total_df_test.mz == total_df_test.mz.max()].mz
			
			max_mass_value = max_mass.to_numpy()

			max_rt = total_df_test[total_df_test.rt == total_df_test.rt.max()].rt
			
			max_rt_value = max_rt.to_numpy()
			self.mMaxBox.setText(str(round(max_mass_value[0])+1))


			self.rtMaxBox.setText(str(round(max_rt_value[0])+1))
			self.mMaxBox.setText(str(round(max_mass_value[0])+1))
	def launch_analysis(self):


		self.window = QtWidgets.QMainWindow()
		self.ui = MainWindowUIClass()
		self.ui.setupUi(self.window)
		self.window.show()
		# self.worker = WorkerThread()
		# self.worker.start()



	def read_me(self):

		os.startfile('README.md')

	def ms1_DF_Export(self):


		if open_state == 2:


			with open((os.path.join(output_path, experiment_name + "_experiment_import_parameters.pickle")), 'rb') as f:
				experiment_info = pickle.load(f)

			if experiment_info.ms2_type == 'DIA':

				print(self.comboBox.currentText())
				comboText=self.comboBox.currentText()
				print(comboText)

				df1, df2, df3 = import_ms1_dataframe(output_path, input_type, comboText)

				df1 = df1[df1.relative_intensity != 0]



				ms2 = import_ms2_dataframe(output_path, input_type, comboText)

				ms2 = ms2[ms2.relative_intensity != 0]




				ms2_data = []
				ms1_data = []

				analyte_list = []

				count = 0
				for analyte_id in df1.replicate_analyte_id.unique():
					ms1_data.append([df1[df1["replicate_analyte_id"] == analyte_id]["relative_intensity"].to_numpy(),
									 df1[df1["replicate_analyte_id"] == analyte_id]["average_mass"].to_numpy()])
					count +=1

				for analyte_id in ms2.replicate_analyte_id.unique():
					ms2_data.append([df1[df1["replicate_analyte_id"] == analyte_id]["relative_intensity"].to_numpy(),
									 df1[df1["replicate_analyte_id"] == analyte_id]["average_mass"].to_numpy()])
					analyte_list.append(analyte_id)
					count +=1



				max_intensity = []
				for analyte in ms1_data:
					res = [x for x in range(len(analyte[0])) if analyte[0][x] == max(analyte[0])] 
					max_intensity.append(max(analyte[1][res]))


				
				filename, _ = QFileDialog.getSaveFileName(self, "Save", os.getcwd()+os.sep+ 'ms_data', "ZIP Files (*.zip)")


				if filename:
					zipObj = ZipFile(filename, 'w')
					i = 0
					for analyte in ms2_data:

						spectrum = Spectrum(mz=analyte[1],
											intensities=analyte[0],
											metadata={"Analyte": analyte_list[i],
													"charge": +1,
													"precursor_mz": max_intensity[i]})

						

						save_as_mgf(spectrum,'analyte_'+ str(analyte_list[i]) + '.mgf')
						zipObj.write('analyte_'+ str(analyte_list[i]) + '.mgf')
						os.remove('analyte_'+ str(analyte_list[i]) + '.mgf')
						i+=1
					
					# if not filename: return 0
					# if filename:


					zipObj.close()

				
			if experiment_info.ms2_type == 'DDA':
				print('test1 ms1 export')
				print(self.comboBox.currentText())
				comboText=self.comboBox.currentText()
				print(comboText)
				print('test2')
				df1, df2, df3 = import_ms1_dataframe(output_path, input_type, comboText)
				print('test3')
				df1 = df1[df1.relative_intensity != 0]
				print('test4')


				ms2 = import_ms2_dataframe(output_path, input_type, comboText)
				print('test5')
				ms2 = ms2[ms2.Intensity != 0]
				print('test6')


				print(df1)
				print(ms2)
				ms2_data = []
				ms1_data = []
				ms1_average_mass = []
				ms2_average_mass = []
				analyte_list = []
				
				count = 0
				for analyte_id in df1.replicate_analyte_id.unique():
					ms1_data.append([df1[df1["replicate_analyte_id"] == analyte_id]["relative_intensity"].to_numpy(),
									 df1[df1["replicate_analyte_id"] == analyte_id]["average_mass"].to_numpy()])
					count +=1
				print(ms1_data)
				for analyte_id in ms2.analyte_id.unique():
					ms2_data.append([ms2[ms2["analyte_id"] == analyte_id]["analyte_id"].to_numpy(),
									 ms2[ms2["analyte_id"] == analyte_id]["Intensity"].to_numpy(),
									 ms2[ms2["analyte_id"] == analyte_id]["ms2_data"].to_numpy(),
									 ms2[ms2["analyte_id"] == analyte_id]["ms1_average_mass"].to_numpy()])
					analyte_list.append(analyte_id)
					ms1_average_mass.append([ms2[ms2["analyte_id"] == analyte_id]["ms1_average_mass"].to_numpy()])
					ms2_average_mass.append([ms2[ms2["analyte_id"] == analyte_id]["ms2_data"].to_numpy()])
					count +=1


				# for mass in ms2.ms1_average_mass.unique():
				# 	ms2_data.append([ms2[ms2["ms1_average_mass"] == mass]["analyte_id"].to_numpy(),
				# 					 ms2[ms2["ms1_average_mass"] == mass]["Intensity"].to_numpy(),
				# 					 ms2[ms2["ms1_average_mass"] == mass]["ms2_data"].to_numpy()])
				# 	ms1_mass_list.append(mass)
				# 	ms1_average_mass.append([ms2[ms2["analyte_id"] == analyte_id]["ms1_average_mass"].to_numpy()])
				# 	ms2_average_mass.append([ms2[ms2["analyte_id"] == analyte_id]["ms2_data"].to_numpy()])
				# 	count +=1


				max_intensity = []
				for analyte in ms1_data:
					res = [x for x in range(len(analyte[0])) if analyte[0][x] == max(analyte[0])] 
					max_intensity.append(max(analyte[1][res]))


				
				filename, _ = QFileDialog.getSaveFileName(self, "Save", os.getcwd()+os.sep+ 'ms_data', "ZIP Files (*.zip)")


				if filename:
					zipObj = ZipFile(filename, 'w')
					i = 0
					for analyte in ms2_data:

						spectrums  = []


						ms2_data_2 = []
						ms1_mass_list = []
						ms2_df_2 = ms2[ms2.analyte_id == analyte_list[i]]
						ms2_df_2 = ms2_df_2.sort_values(by=['ms2_data'])
						for mass in ms2_df_2.ms1_average_mass.unique():
							ms2_data_2.append([ms2_df_2[ms2_df_2["ms1_average_mass"] == mass]["ms2_data"].to_numpy(),
											   ms2_df_2[ms2_df_2["ms1_average_mass"] == mass]["Intensity"].to_numpy()])

							ms1_mass_list.append(mass)

						print('ms1_mass_list', ms1_mass_list)
						print('ms2_data_2', ms2_data_2)
						print('MS2 DF 2',ms2_df_2)
						j=0
						for mass in ms2_data_2:
							print('i',i)
							print('j',j)
							print('analyte', analyte_list[i])
							print('ms1_mass', ms1_mass_list[j])
							print(mass)
							print(mass[0])
							ms2_mass = mass[0]
							ms2_intensity = mass[1]
							print(ms2_mass)
							print(mass[1])

							# spectrum = Spectrum(mz=numpy.array(mass[0],dtype="float"),
							# 					intensities=numpy.array(mass[1],dtype="float"),
							# 					metadata={"Analyte": analyte_list[i],
							# 							"charge": +1,
							# 							"precursor_mz": ms1_mass_list[j]})
							spectrum = Spectrum(mz=numpy.array(ms2_mass, dtype="float"),
												intensities=numpy.array(ms2_intensity, dtype="float"),
												metadata={"Analyte": analyte_list[i],
														"charge": +1,
														"precursor_mz": ms1_mass_list[j]})

							print('test')

							spectrums.append(spectrum)
							j+=1
						save_as_mgf(spectrums,'analyte_'+ str(analyte_list[i]) + '.mgf')
						zipObj.write('analyte_'+ str(analyte_list[i]) + '.mgf')
						os.remove('analyte_'+ str(analyte_list[i]) + '.mgf')
							
						i+=1
					
					# if not filename: return 0
					# if filename:


					zipObj.close()




		return


	def ms2_Plot_DF(self):


		print('ms2_plot test 1')

		if open_state == 2:
			print('ms2_plot test 1')
			comboText=self.comboBox.currentText()
			print('ms2_plot test 3')
			with open((os.path.join(output_path, experiment_name + "_experiment_import_parameters.pickle")), 'rb') as f:
				experiment_info = pickle.load(f)
			print('ms2_plot test 4')





			if experiment_info.ms2_type == 'DIA':
				print('ms2_plot test 5')
				ms2 = import_ms2_dataframe(input_structure, input_type, comboText)

				print('ms2_plot test 6')
				print('MS2 DF',ms2)

				ms2.sort_values('average_mass', inplace=True)
				print('ms2_plot test 7')



			if experiment_info.ms2_type == 'DDA':
				ms2 = import_ms2_dataframe(input_structure, input_type, comboText)


				ms2.sort_values('ms2_data', inplace=True)

				print('==========================================================================',ms2[ms2.analyte_id == 7])
		else:
			ms2 = []
		print('ms2', ms2)
		return ms2

	def nullButton(self):
	    # method called by button 

		# if button is checked 
		if self.Nullbutton.isChecked():
			#self.update_tab()
			# setting background color to light-blue 
			self.Nullbutton.setStyleSheet("background-color : lightgrey") 
			self.Nullbutton.setText('Off') 

		# if it is unchecked 
		else: 
			#self.update_tab()
			# set background color back to light-grey 
			self.Nullbutton.setStyleSheet("background-color : lightgreen") 
			self.Nullbutton.setText('On')


		
		



	def blankButton(self):

	    # method called by button 

	  
	    # if button is checked 
		if self.Blankbutton.isChecked():
			#self.update_tab()
	        # setting background color to light-blue 
			self.Blankbutton.setStyleSheet("background-color : lightgrey") 
			self.Blankbutton.setText('Off') 

	    # if it is unchecked 
		else: 
			#self.update_tab()
			# set background color back to light-grey 
			self.Blankbutton.setStyleSheet("background-color : lightgreen") 
			self.Blankbutton.setText('On')
	def toggle_rt(self):

	    # method called by button 

	  
	    # if button is checked 
		if self.Toggle_rt.isChecked():
			#self.update_tab()
	        # setting background color to light-blue 
			self.Toggle_rt.setStyleSheet("background-color : lightgrey")
			self.Toggle_rt.setText('Toggle Line Plot') 

	    # if it is unchecked 
		else: 
			#self.update_tab()
			# set background color back to light-grey 
			self.Toggle_rt.setStyleSheet("background-color : lightgrey")
			self.Toggle_rt.setText('Toggle Scatter Plot')

	def mMinLimit(self, text):

		print('minlimit',text)




		maxvalue=float(self.mMaxBox.text())
		if ',' in text:
			text = 0.0
			self.mMinBox.setStyleSheet("background:rgb(255,0,0);")
		else:
			self.mMinBox.setStyleSheet("background:rgb(255,255,255);")
		if text == "":
			text = 0.0
			self.mMinBox.setText('0')
			if maxvalue == "":
				maxvalue = 10.0

		if maxvalue == "":
			maxvalue = 10.0
			if text == "":
				text = 0.0
				self.mMinBox.setText('0')

		print('mMinLimit',text)
		print('MaxValue',maxvalue)
		self.mz_plot.setXRange(min=float(text),max=float(maxvalue))
		self.mz_r1_plot.setXRange(min=float(text),max=float(maxvalue))
		self.mz_r2_plot.setXRange(min=float(text),max=float(maxvalue))
		self.mz_r3_plot.setXRange(min=float(text),max=float(maxvalue))

	def mMaxLimit(self, text):

		print('maxlimit',text)


		minvalue=float(self.mMinBox.text())
		if ',' in text:
			text = 0.0
			self.mMaxBox.setStyleSheet("background:rgb(255,0,0);")
		else:
			self.mMaxBox.setStyleSheet("background:rgb(255,255,255);")

		print(minvalue)

		if text == "":
			text = 0
			self.mMaxBox.setText('0')
			if minvalue == "":
				minvalue = 0.0


		if minvalue == "":
			minvalue = 0.0
			if text == "":
				text = 0
				self.mMaxBox.setText('0')


		print('mMaxLimit',text)
		print('Min Value',minvalue)
		self.mz_plot.setXRange(min=float(minvalue),max=float(text))
		self.mz_r1_plot.setXRange(min=float(minvalue),max=float(text))
		self.mz_r2_plot.setXRange(min=float(minvalue),max=float(text))
		self.mz_r3_plot.setXRange(min=float(minvalue),max=float(text))




	def rtMinLimit(self, text):

		print('minlimit',text)




		maxvalue=float(self.rtMaxBox.text())
		if ',' in text:
			text = 0.0
			self.rtMinBox.setStyleSheet("background:rgb(255,0,0);")
		else:
			self.rtMinBox.setStyleSheet("background:rgb(255,255,255);")
		if text == "":
			text = 0.0
		if text == "":
			text = 0.0
			self.rtMinBox.setText('0')
			if maxvalue == "":
				maxvalue = 10.0

		if maxvalue == "":
			maxvalue = 10.0
			if text == "":
				text = 0.0
				self.rtMinBox.setText('0')

		print('mMinLimit',text)
		print('MaxValue',maxvalue)




		# if text == "":
		# 	text = 0

 
		# print(self.rtMaxBox.text())
		# maxvalue=float(self.rtMaxBox.text())
		# if maxvalue == "":
		# 	maxvalue = 0
		self.rt_plot.setXRange(min=float(text),max=maxvalue)
		self.rt_r1_plot.setXRange(min=float(text),max=maxvalue)
		self.rt_r2_plot.setXRange(min=float(text),max=maxvalue)
		self.rt_r3_plot.setXRange(min=float(text),max=maxvalue)


	def rtMaxLimit(self, text):

		print('maxlimit',text)


		minvalue=float(self.rtMinBox.text())
		if ',' in text:

			text = 0.0
			self.rtMaxBox.setStyleSheet("background:rgb(255,0,0);")
		else:
			self.rtMaxBox.setStyleSheet("background:rgb(255,255,255);")

		print(minvalue)

		if text == "":
			text = 10.0
			self.rtMaxBox.setText('0')
			if minvalue == "":
				minvalue = 0.0


		if minvalue == "":
			minvalue = 0.0
			if text == "":
				text = 10.0
				self.rtMaxBox.setText('0')


		print('mMaxLimit',text)
		print('Min Value',minvalue)


		# if text == "":
		# 	text = 0
		# print(self.rtMinBox.text())
		# minvalue=float(self.rtMinBox.text())

		# if minvalue == "":
		# 	minvalue = 0
		self.rt_plot.setXRange(min=minvalue,max=float(text))
		self.rt_r1_plot.setXRange(min=minvalue,max=float(text))
		self.rt_r2_plot.setXRange(min=minvalue,max=float(text))
		self.rt_r3_plot.setXRange(min=minvalue,max=float(text))


	def hide_sample(self):

		#self.btn2.setEnabled(True)
		self.mz_plot.hide()
		self.rt_plot.hide()
		#self.ms1_plot.hide()
		#self.ms2_plot.hide()

	def show_sample(self):

		tab_state = Tab_State(self.btn2.setEnabled(False), self.btn3.setEnabled(True), self.btn4.setEnabled(True), self.btn5.setEnabled(True))
		self.btn2.setEnabled(False)
		self.btn4.setEnabled(True)
		self.btn5.setEnabled(True)
		self.btn3.setEnabled(True)
		# self.Nullbutton.setEnabled(True)
		# self.chooseSample.setEnabled(True)
		# self.showNull.setEnabled(True)
		#self.comboBox.setEnabled(True)
		self.ms1_plot.show()
		self.ms2_plot.show()


		self.mz_plot.show()
		self.rt_plot.show()
			#self.ms1_plot.show()
			#self.ms2_plot.show()
		self.btn2.setStyleSheet('background-color : lightgreen')
		self.btn5.setStyleSheet('background-color : lightgrey')
		self.btn3.setStyleSheet('background-color : lightgrey')
		self.btn4.setStyleSheet('background-color : lightgrey')
		
		
	def show_replicate(self):
		tab_state = Tab_State(self.btn2.setEnabled(True), self.btn3.setEnabled(False), self.btn4.setEnabled(True), self.btn5.setEnabled(True))
		self.btn3.setEnabled(False)

		self.btn2.setEnabled(True)
		self.btn5.setEnabled(True)
		self.btn4.setEnabled(True)
		# self.Nullbutton.setEnabled(True)
		# self.chooseSample.setEnabled(True)
		# self.showNull.setEnabled(True)

		self.mz_r1_plot.show()
		self.mz_r2_plot.show()
		self.mz_r3_plot.show()
		self.rt_r1_plot.show()
		self.rt_r2_plot.show()
		self.rt_r3_plot.show()
		self.ms1_plot.show()
		self.ms2_plot.show()

		self.btn3.setStyleSheet('background-color : lightgreen')
		self.btn5.setStyleSheet('background-color : lightgrey')
		self.btn2.setStyleSheet('background-color : lightgrey')
		self.btn4.setStyleSheet('background-color : lightgrey')


	def hide_replicate(self):
		self.mz_r1_plot.hide()
		self.mz_r2_plot.hide()
		self.mz_r3_plot.hide()
		self.rt_r1_plot.hide()
		self.rt_r2_plot.hide()
		self.rt_r3_plot.hide()

	def show_experiment(self):

		print('show experiment test 1')
		tab_state = Tab_State(self.btn2.setEnabled(True), self.btn3.setEnabled(True), self.btn4.setEnabled(False), self.btn5.setEnabled(True))
		self.btn4.setEnabled(False)


		self.btn2.setEnabled(True)
		self.btn3.setEnabled(True)
		self.btn5.setEnabled(True)
		# self.Nullbutton.setEnabled(False)
		self.ex_bOn_plot.show()
		# self.chooseSample.setEnabled(False)
		# self.showNull.setEnabled(False)
		self.ms1_plot.show()
		self.ms2_plot.show()
		# if self.Blankbutton.isChecked() == False:
		# 	self.ex_bOn_plot.show()
		# 	#self.ex_bOff_plot.hide()
		# if self.Blankbutton.isChecked()  == True:
		# 	self.ex_bOn_plot.hide()
			#self.ex_bOff_plot.show()
		self.btn4.setStyleSheet('background-color : lightgreen')

		self.btn5.setStyleSheet('background-color : lightgrey')
		self.btn3.setStyleSheet('background-color : lightgrey')
		self.btn2.setStyleSheet('background-color : lightgrey')

	def hide_experiment(self):

		self.ex_bOn_plot.hide()
		#self.ex_bOff_plot.hide()

	def show_diversity(self):
		tab_state = Tab_State(self.btn2.setEnabled(True), self.btn3.setEnabled(True), self.btn4.setEnabled(True), self.btn5.setEnabled(False))
		self.btn4.setEnabled(True)
		self.btn2.setEnabled(True)
		self.btn3.setEnabled(True)
		self.btn5.setEnabled(False)
		self.btn5.setStyleSheet('background-color : lightgreen')

		self.btn4.setStyleSheet('background-color : lightgrey')
		self.btn3.setStyleSheet('background-color : lightgrey')
		self.btn2.setStyleSheet('background-color : lightgrey')
		# self.Nullbutton.setEnabled(True)
		# self.chooseSample.setEnabled(True)
		# self.showNull.setEnabled(True)
		self.div_plot.show()
	def hide_diversity(self):
		self.div_plot.hide()
		pass







	def close_application(self):
		choice = QMessageBox.question(self, 'Quit',
											"Are you sure you want to quit?",
											QMessageBox.Yes | QMessageBox.No)
		if choice == QMessageBox.Yes:
			print("Closing Application")
			sys.exit()
		else:
			pass
			
# class WorkerThread(QThread):
# 	def run(self):
# 		self.window = QtWidgets.QMainWindow()
# 		self.ui = MainWindowUIClass()
# 		self.ui.setupUi(self.window)
# 		self.window.show()

		
def run():

	app = QApplication(sys.argv)
	app.setStyle('Fusion')
	GUI = Window()
	sys.exit(app.exec_())
run()
