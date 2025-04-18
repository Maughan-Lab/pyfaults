# initial gui for faults

# imports
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import sys
import os
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import time

# Import pyfaults!
import pyfaults
from pyfaults import * 
from pyfaults.pfInput import pfInput  # If pfInput is a function or class inside the pfInput.py file
from pyfaults.simXRD import fullSim
from pyfaults import tt_to_q


# =========================================================================================================================================
# MAIN WINDOW
# =========================================================================================================================================
class MainWindow(QMainWindow):
    # =====================================================================================================================================
    #  INSTANCE VARIABLES
    # =====================================================================================================================================

    wavelength_inputs = []
    
    # TODO: refactor to include "file_path" or similar to name
    diff_data = "" #contains file path to diffraction data
    obs_data = "" #contains file path to observed data
    diff_cif_fp = "" #contains file path to csv
    
    status_msg = "" #contains message about status that appears beside Generate button
    # TODO: remove global variable (not needed)
    
    # create a dictionary to hold all the necessary fault information
    fault_info = {
        "Layer": 0, 
        "N": 0,
        "2Theta": 0,
        "Wavelength": 0,
        "Broadening": 0,
        "Vector": [0,0,0],
        "Probability": []
    }
    # create a pandas DF to contain xy data from Diffraction Pattern and Observed Data
    diffx_df = pd.DataFrame(columns=["x", "y"])
    obs_df = pd.DataFrame(columns=["x", "y"])
    
    # create a variable to check if the graph has been generated
    is_graphed = False

    # ===================================================================================================================================
    # INITIALIZATION
    # ===================================================================================================================================
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setWindowTitle("Faults Form") #window title
        self.setGeometry(0,0,1000,800) #this shrinks the size when clicked on the minimize
        # set style of the window
        self.setStyleSheet(
            "background: 'white';" +
            "font-family: 'helvetica'; " 
        )

        # creating a central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        # creating a grid layout
        self.grid = QGridLayout() #currently row: 15 (TODO: adjust if needed), col: 20
        # set stretch factors to the grid to dynamically adjust row and col spans to fit screen size
        # number of rows will need to be adjusted after determining how many inputs will go into the visualizer
        for row in range(15):
            self.grid.setRowStretch(row, 1)
        for col in range(20):
            self.grid.setColumnStretch(col, 1)
        
        
        # create labels for the top of the page
        # left page label - GUI title
        self.title = QLabel("Pyfaults Visualizer")
        self.title.setAlignment(Qt.AlignLeft)
        # right page label - name of graphs (TODO: update when clicking generate button)
        self.graph_title = QLabel("Graph Title")
        self.graph_title.setAlignment(Qt.AlignLeft)
        
        # add line to split
        self.split_line = QFrame()
        self.split_line.setFrameShape(QFrame.VLine) #vertical line
        self.split_line.setStyleSheet(
            "border: 2px solid #eeeeee;"
        )     
                    
        # add both labels to the grid
        self.grid.addWidget(self.title, 0, 0, 1, 8)
        self.grid.addWidget(self.split_line, 0, 8, 15, 1) #adding split line in middle as well
        self.grid.addWidget(self.graph_title, 0, 9, 1, 11)
 
 
        # -- UPLOAD BUTTONS --
        upload_icon = QIcon(QPixmap("visualizer_gui/icons/upload.png").scaled(40, 40)) #TODO: make scaled dyanmically
        
        # start with upload files - both are buttons that will call another function that will allow files to be uploaded
        # diffraction pattern upload
        self.upload_diff = QPushButton("Diffraction Pattern")
        self.upload_diff.setStyleSheet(
            "*{" +
            "text-align: left;" +
            "color: '#333333';" +
            "font-family: 'helvetica';" +
            "border: 0.05em solid 'gray';" +
            "border-radius: 0.5em;" +
            "padding: 0.5em;" +
            "margin: 0.5em;" +
            "}*" +
            ":hover{background: '#EEEEEE';}" + #background when hovering over button
            ":pressed{background: '#DDDDDD';}" #even darker background when button clicked
        )
        self.upload_diff.clicked.connect(lambda: self.upload_file("Text Files (*.txt)")) #connect button to upload_file() method
        self.upload_diff.setIcon(upload_icon) #include upload icon


        # measured values file upload
        self.upload_obs = QPushButton("Observed Data")
        self.upload_obs.setStyleSheet(
            "*{" +
            "text-align: left;" +
            "color: '#333333';" +
            "font-family: 'helvetica';" +
            "border: 0.05em solid 'gray';" +
            "border-radius: 0.5em;" +
            "padding: 0.5em;" +
            "margin: 0.5em;" +
            "}*" +
            ":hover{background: '#EEEEEE';}" + #background when hovering over button
            ":pressed{background: '#DDDDDD';}" #even darker background when button clicked
        )
        self.upload_obs.clicked.connect(lambda: self.upload_file("XY Files (*.xy)")) #connect button to upload_file() method
        self.upload_obs.setIcon(upload_icon) #include upload icon
       
        # add the upload buttons to the grid
        self.grid.addWidget(self.upload_diff, 1, 0, 1, 4)
        self.grid.addWidget(self.upload_obs, 1, 4, 1, 4)
             
        
        # -- FAULTS INFO INPUT
        # TODO: style the input spaces with style sheet
        self.faults_info_header = QLabel("Faults Information")
        # layer
        self.layer_info = QLineEdit()
        self.layer_info.setPlaceholderText("Layer")
        # N
        self.n_info = QLineEdit()
        self.n_info.setPlaceholderText("N")
        # 2Theta
        self.theta_info = QLineEdit()
        self.theta_info.setPlaceholderText("2\u03B8")
        # Wavelength
        self.wavelength_info = QLineEdit()
        self.wavelength_info.setPlaceholderText("Wavelength")
        # Broadening
        self.broadening_info = QLineEdit()
        self.broadening_info.setPlaceholderText("Broadening")
        # vector
        self.vector_header = QLabel("Vector")
        self.vector_x = QLineEdit()
        self.vector_x.setPlaceholderText("x")
        self.vector_y = QLineEdit()
        self.vector_y.setPlaceholderText("y")
        self.vector_z = QLineEdit()
        self.vector_z.setPlaceholderText("z")
        # probability
        self.probability_header = QLabel("Probability")
        self.prob_info = QLineEdit()
        self.prob_info.setPlaceholderText("p1, p2, ...")
        
        # change x min and x max 
        self.graph_info = QLabel("Graph Information")
        self.x_min = QLineEdit()
        self.x_min.setPlaceholderText("X Min (Optional)")
        self.x_max = QLineEdit()
        self.x_max.setPlaceholderText("X Max (Optional)")
        
        self.diffx_color_info = QLabel("Diffraction Color:")
        self.diffx_color = QComboBox()
        self.diffx_color.addItems(["Red", "Blue", "Green", "Orange", "Black"])
        self.diffx_color.setCurrentText("Red")  # Default to Red
        self.obs_color_info = QLabel("Observed Color:")
        self.obs_color = QComboBox()
        self.obs_color.addItems(["Red", "Blue", "Green", "Orange", "Black"])
        self.obs_color.setCurrentText("Blue")  # Default to Blue
        self.graph_update_button = QPushButton("Update Graph")
    
        
        # add all input information into grid
        self.grid.addWidget(self.faults_info_header, 2, 0, 1, 4)
        self.grid.addWidget(self.layer_info, 3, 0, 1, 1)
        self.grid.addWidget(self.n_info, 3, 1, 1, 1)
        self.grid.addWidget(self.theta_info, 3, 2, 1, 1)
        self.grid.addWidget(self.wavelength_info, 3, 3, 1, 2)
        self.grid.addWidget(self.broadening_info, 3, 5, 1, 2)
        self.grid.addWidget(self.vector_header, 4, 0, 1, 4)
        self.grid.addWidget(self.vector_x, 5, 0, 1, 1)
        self.grid.addWidget(self.vector_y, 5, 1, 1, 1)
        self.grid.addWidget(self.vector_z, 5, 2, 1, 1)
        self.grid.addWidget(self.probability_header, 4, 3, 1, 4)
        self.grid.addWidget(self.prob_info, 5, 3, 1, 4)
        self.grid.addWidget(self.graph_info, 6, 0, 1, 4)
        self.grid.addWidget(self.x_min, 7, 0, 1, 2)
        self.grid.addWidget(self.x_max, 7, 2, 1, 2)
        self.grid.addWidget(self.diffx_color_info, 8, 0, 1, 2)
        self.grid.addWidget(self.diffx_color, 8, 2, 1, 2)
        self.grid.addWidget(self.obs_color_info, 9, 0, 1, 2)
        self.grid.addWidget(self.obs_color, 9, 2, 1, 2)
        
        # -- BOTTOM BUTTONS --
        
        # update graph button
        self.graph_update_button.setStyleSheet(
            "*{" +
            "text-align: center;" +
            "color: 'white';" +
            "background: 'black';" +
            "font-family: 'helvetica';" +
            "border: 0.05em solid 'black';" +
            "border-radius: 0.5em;" +
            "padding: 0.5em;" +
            "margin: 0.5em;" +
            "}*" +
            ":hover{background: '#222222';}" + #background when hovering over button
            ":pressed{background: '#333333';}" #even darker background when button clicked
        )
        self.graph_update_button.clicked.connect(self.update_graphs)
        self.grid.addWidget(self.graph_update_button, 6, 6, 1, 2) #add button to grid
        
        # generate button
        self.generate_button = QPushButton("Generate")
        self.generate_button.setStyleSheet(
            "*{" +
            "text-align: center;" +
            "color: 'white';" +
            "background: 'black';" +
            "font-family: 'helvetica';" +
            "border: 0.05em solid 'black';" +
            "border-radius: 0.5em;" +
            "padding: 0.5em;" +
            "margin: 0.5em;" +
            "}*" +
            ":hover{background: '#222222';}" + #background when hovering over button
            ":pressed{background: '#333333';}" #even darker background when button clicked
        )
        self.generate_button.clicked.connect(self.generate) #connect button to generate() method 
        # TODO: add restriction where this can't occur until file is uploaded
        self.grid.addWidget(self.generate_button, 14, 0, 1, 2) #add button to grid
     
    
        # export button 
        self.export_button = QPushButton("Export")
        self.export_button.setStyleSheet(
            "*{" +
            "text-align: center;" +
            "color: 'white';" +
            "background: 'black';" +
            "font-family: 'helvetica';" +
            "border: 0.05em solid 'black';" +
            "border-radius: 0.5em;" +
            "padding: 0.5em;" +
            "margin: 0.5em;" +
            "}*" +
            ":hover{background: '#222222';}" + #background when hovering over button
            ":pressed{background: '#333333';}" #even darker background when button clicked
        )
        # TODO: connect to save function to save both the .xy file from diffraction pattern and save the graphs
        self.grid.addWidget(self.export_button, 14, 6, 1, 2) #add button to grid
        
        
        # -- STATUS MESSAGE --
        # status message - add in a space to place message about status
        self.status = QLabel(self.status_msg)
        self.status.setStyleSheet(
            "text-align: left;" +
            "color: 'black';" +
            "font-family: 'helvetica';" 
        )
        self.status.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)  # Let label expand horizontally
        # To ensure that text truncates with "..." when needed
        self.status.setWordWrap(True)  # Disable word wrapping
        self.grid.addWidget(self.status, 14, 2, 1, 4) #add message to grid - initialized to empty
        
        
        # add set the window layout to the grid
        central_widget.setLayout(self.grid)
        
        
        # TODO: create a help menu and help button with shortcuts (Ctrl+?)
               
        # Ctrl+W close application
        close_shortcut = QShortcut(QKeySequence("Ctrl+W"), self)
        close_shortcut.activated.connect(self.close_application)




    # =========================================================================================================================================
    # FUNCTIONS
    # =========================================================================================================================================

    def close_application(self):
        print("Closing Application") #TODO: have this check to save files before closing
        self.close()


    # Generating Correct Files and Graphs
    def generate(self):
        # first check if both files are uploaded
        if not self.diff_data:  # Ensure the file path is set for diffractin patter data
            self.status_msg = "ERROR: No file selected for Diffraction Pattern."
            self.status.setText(self.status_msg)
            self.status.setStyleSheet("color: red; font-weight: bold;") #set error message to red and bold
            return
        if not self.obs_data:  # Ensure the file path is set for observed data
            self.status_msg = "ERROR: No file selected for Observed Data."
            self.status.setText(self.status_msg)
            self.status.setStyleSheet("color: red; font-weight: bold;") #set error message to red and bold
            return
        # update fault info with values inputted
        self.update_fault_info()
               
        # check if fault info is valid - calls the is_fault_info_valid() method
        if not self.is_fault_info_valid():
            self.status_msg = "ERROR: Invalid input in Faults Information section."
            self.status.setText(self.status_msg)
            self.status.setStyleSheet("color: red; font-weight: bold;") #set error message to red and bold
            return 

        
        # after checking necessary values added -> use pyfaults to read diffraction pattern
        start_time = datetime.now()
        self.status_msg = "Reading Diffraction Patterns"
        self.status.setText(self.status_msg)
        self.status.setStyleSheet("color: black; font-weight: normal;") 
        # QTimer.singleShot(1000, self.update_status("Reading Diffraction Patters"))
        # call method to read through Diffraction Pattern file
        self.read_diff_data()

        # simulating supercell with faults info
        self.status_msg = "Simulating Supercell"
        self.status.setText(self.status_msg)
        self.status.setStyleSheet("color: black; font-weight: normal;") 
        # call method to simulate 
        self.simulate_fault()
        
        
        # normalizing observed data
        self.status_msg = "Normalizing Data"
        self.status.setText(self.status_msg)
        self.status.setStyleSheet("color: black; font-weight: normal;") 
        # call method to normalize obs data using tt_to_q from pyfaults
        self.normalize_obs_data()
        self.normalize_diff_data()
        
        # creating graphs
        self.status_msg = "Creating Graphs"
        self.status.setText(self.status_msg)
        self.status.setStyleSheet("color: black; font-weight: normal;") 
        # create graph method 
        self.create_graph()
        
        # after creating plots - send message of completing process and then output time
        end_time = datetime.now()
        # total_seconds = (end_time - start_time).total_seconds()
        # hrs = int(total_seconds // 3600)  # Number of full hours
        # min = int((total_seconds % 3600) // 60)  # Remaining minutes
        # sec = int(total_seconds % 60)  # Remaining seconds
        # self.status_msg = f"COMPLETE: Time Elapsed {hrs}hrs {min}min {sec}s"
        self.status_msg = f"COMPLETE: Time Completed - {end_time}"
        self.status.setText(self.status_msg)
        self.status.setStyleSheet("color: green; font-weight: bold;") 
        
    # function to update status message
    # TODO: fix delay to work properly
    def update_status(self, msg):
        self.status.setText(msg)
        self.status.setStyleSheet("color: black; font-weight: normal;") 
        
    # function to update fault info values 
    def update_fault_info(self):
        # TODO: check that all the values are filled - this might be caught in the "is_fault_info_valid()" function instead (ex. data type match, valid range, etc.)
        
        # check that all fields are not empty before trying to update fields
        if (self.layer_info.text() and self.n_info.text() and self.theta_info.text() and self.wavelength_info.text() and self.broadening_info.text() and 
            self.vector_x.text() and self.vector_y.text() and self.vector_z.text() and self.prob_info.text()):
            try:
                new_values = {
                    "Layer": int(self.layer_info.text()),
                    "N": int(self.n_info.text()),
                    "2Theta": float(self.theta_info.text()),
                    "Wavelength": float(self.wavelength_info.text()),
                    "Broadening": float(self.broadening_info.text()),
                    "Vector": [float(self.vector_x.text()), float(self.vector_y.text()), float(self.vector_z.text())],
                    "Probability": [float(x) for x in self.prob_info.text().split(",")]
                }
                self.fault_info.update(new_values)
            except ValueError as e:
                # Handle the error if casting fails (e.g., non-numeric input)
                self.status_msg = f"ERROR: Non-numeric value detected in the Faults Information section."
                self.status.setText(self.status_msg)
                self.status.setStyleSheet("color: red; font-weight: bold;")
        
    # function to check if all the values inputted are valid - returns boolean
    def is_fault_info_valid(self):
        if (self.fault_info["Layer"] < 1) : return False
        if (self.fault_info["N"] < 1) : return False
        if (len(self.fault_info["Probability"]) < 1): return False
        # each of the probabilities should be between 0-1
        # each of vector x,y,z should be between 0-1(not including 0 and 1)
        #2 theta max should be 120
        # broadening should also be positive
        # TODO: add more conditions
        return True
        
        
    # read diffracion data file
    def read_diff_data(self):  
        # TODO: make sure in every diffraction data file to remove the word "Displacement" from the "TYPE:" field to prevent calling the genSupercells() method
        # For now, use the "edited" version of the txt file

        # create unit cells
        unitcell, ucDF, gsDF, scDF, simDF = pfInput(self.diff_data) #using function in pfinput.py
        #unitcell = unitcell object
        #ucDF = df containing information about layers 

        # create super cells from the unit cells
        # print(self.fault_info) #test print
        final_supercell = Supercell(unitcell, self.fault_info["N"], "N"+ str(self.fault_info["Layer"]), self.fault_info["Vector"], self.fault_info["Probability"][0])
        '''
        Parameters
        ----------
        unitcell (Unitcell) : base unit cell
        nStacks (int) : number of unit cell stacks in supercell
        fltLayer (str, optional) : faulted layer name
        stackVec (array_like, optional) : displacement vector [x,y,z] for faulted layer in fractional coordinates
        stackProb (float, optional) : stacking fault probability
        '''
        
        # turn into cif file and store the path in the corresponding variable
        # get directory path
        dir_path = os.path.dirname(self.diff_data) 
        # Add a trailing slash if it doesn't already have one
        if not dir_path.endswith(os.sep):
            dir_path += "/"
        # test print directory path
        # print("dir path", dir_path)
        # turn to cif using .toCif() method in Supercell class
        final_supercell.toCif(dir_path) #pass dir 
        
        # save the cif dir path to global var
        # NOTE: if the .toCif() method changes in the Supercell class, this following code might have to change as well
        self.diff_cif_fp = dir_path + 'Supercell_' + unitcell.name + '_N' + str(self.fault_info["N"]) + '_CIF.cif'
        
        
        
        # TODO: Implement the following - completed above -> remove comments when appropriate
        # remove the simulation type (ex. Displacement) from LYC_input.txt file so that it will take in the info from the txt file
        # take this information from unit cell and potentially turn it into a cif file?
        # take the inputs from the gui and edit the text file fill it with this info
        # pfinput.py -> get unit cell
        # unit cell -> create super cell
        # super cell -> output into cif file
        # The following is done in the simulate_faults() method
        # cif file -> do the diffx 
        # diffx q and ints -> plots
    
    def simulate_fault(self):
        # calling the XRD method for a full simulation 
        dir_path = os.path.dirname(self.diff_cif_fp)
        # Add a trailing slash if it doesn't already have one
        if not dir_path.endswith(os.sep):
            dir_path += "/"
        file_name = os.path.splitext(os.path.basename(self.diff_cif_fp))[0] #gets file name only and removes extension
        print("Directory path:", dir_path)
        print("File name:", file_name)
        
        # run full simulation
        '''
        Parameters
        ----------
        path (str) : CIF file location
        cif (str) : CIF file name
        wl (float) : instrument wavelength
        tt_max (float) : maximum two theta
        savePath (str) : location to save simulation data to
        pw (float, optional) : peak broadening term
        bg (float, optional) : average of normal background

        Returns
        -------
        q (array_like) : calculated Q values
        ints (array_like) : calculated intensity values, normalized
        '''
        # simulates XRD pattern, returns normalized intensity values using fullSim() function in simXRD pyfaults file
        # TODO: change back to dir_path eventually
        fullSim(dir_path, file_name, self.fault_info["Wavelength"], self.fault_info["2Theta"], dir_path, pw=self.fault_info["Broadening"])
        # save values in diffraction data df
        self.sim_fp = dir_path + file_name + '_sim.txt'
        sim_data = pd.read_csv(self.sim_fp, delim_whitespace=True, header=None, names=["x", "y"])
        self.diffx_df["x"] = sim_data["x"]
        self.diffx_df["y"] = sim_data["y"]
            
    # method to normalize observed data
    def normalize_obs_data(self):
        # Load the data from the .xy file (assuming the file is space or tab-separated)
        data = np.loadtxt(self.obs_data)
        
        # Extract x and y values (assuming the file has two columns)
        self.obs_df["x"] = data[:, 0]
        self.obs_df["x"] = tt_to_q(self.obs_df["x"], self.fault_info["Wavelength"])
        self.obs_df["y"] = data[:, 1]
        
        # call function to change x-axis from tt to q
        # self.obs_data["x"] = np.array(pyfaults.tt_to_q(data[:,0], self.fault_info["Wavelength"]))

        # 0 to 1 normalization
        min_val = self.obs_df["y"].min()
        max_val = self.obs_df["y"].max()
        self.obs_df["y"] = (self.obs_df["y"] - min_val) / (max_val - min_val)

    # method to normalize diffraction data (but should this already come out normalized?)
    def normalize_diff_data(self):
        # data contained in self.diffx_df
        
        # 0 to 1 normalization
        min_val = self.diffx_df["y"].min()
        max_val = self.diffx_df["y"].max()
        self.diffx_df["y"] = (self.diffx_df["y"] - min_val) / (max_val - min_val)
    
    # upload file method
    # TODO: update method with params to accept only certain file types
    def upload_file(self, file_types="All Files (*)"):
        '''
        file_types: specify what file type is allowed - defaulted to all types
        '''
        # Open a file dialog and allow the user to select a file
        # Define options for file dialong
        options = QFileDialog.Options() 
        options |= QFileDialog.ReadOnly #Add ReadOnly option for file
        # call function to retrieve file path
        file_path, _ = QFileDialog.getOpenFileName(
            self, 
            "Select a File to Upload", 
            "", 
            file_types, #allowed file types passed as param by the button that called function
            options=options #set options to the read only from above
        )
        # If a file is selected, display its path in label below
        if file_path:
            button = self.sender() #get button that triggered this function
            # set this to appropriate global variable
            if button == self.upload_diff: #if button was "Diffraction Pattern" -> save file path to diff_data 
                self.diff_data = file_path
                self.fill_data()
            elif button == self.upload_obs: #if button was "Observed Data" -> save file path to obs_data 
                self.obs_data = file_path
            # get file name from the path and replace button name
            file_name = os.path.basename(file_path)
            # set text of that button if a file was given
            button.setText(f"File: {file_name}")


    # this function will fill in current fields using the information in the diffraction pattern file
    def fill_data(self):
        # first find all the information in the file and fill in the fault_info dictionary, if you can't find it, leave it empty
        # Open and read the file
        with open(self.diff_data, "r") as file:
            for l in file:
                line = l.strip()
                
                # Extract Layer
                if line.startswith("FAULT LAYER:L"):
                    self.fault_info["Layer"] = int(line[len("FAULT LAYER:L"):].strip())

                # Extract N
                if line.startswith("N:"):
                    self.fault_info["N"] = int(line[2:].strip())

                # Extract Probability (comma-separated values)
                if line.startswith("PROBABILITY:"):
                    prob_values = line[len("PROBABILITY:"):].strip().split(",")
                    self.fault_info["Probability"] = [float(p) for p in prob_values]

                # Extract Vector (comma-separated inside brackets)
                if line.startswith("VECTOR:[") and line.endswith("]"):
                    vector_values = line[8:-1].strip().split(",")
                    self.fault_info["Vector"] = [float(v) for v in vector_values]

                # Extract Wavelength
                if line.startswith("WAVELENGTH:"):
                    self.fault_info["Wavelength"] = float(line[len("WAVELENGTH:"):].strip())

                # Extract 2Theta
                if line.startswith("MAX TWO THETA:"):
                    self.fault_info["2Theta"] = float(line[len("MAX TWO THETA:"):].strip())

                # Extract Broadening
                if line.startswith("BROADENING:"):
                    self.fault_info["Broadening"] = float(line[len("BROADENING:"):].strip())
        
        # second use this dictionary to fill in the text in the GUI fields
        self.layer_info.setText(str(self.fault_info["Layer"])) 
        self.n_info.setText(str(self.fault_info["N"])) 
        self.theta_info.setText(str(self.fault_info["2Theta"]))  
        self.wavelength_info.setText(str(self.fault_info["Wavelength"]))  
        self.broadening_info.setText(str(self.fault_info["Broadening"])) 
        self.vector_x.setText(str(self.fault_info["Vector"][0]))  
        self.vector_y.setText(str(self.fault_info["Vector"][1]))  
        self.vector_z.setText(str(self.fault_info["Vector"][2]))  
        self.prob_info.setText(str(self.fault_info["Probability"][0])) 


    def create_graph(self, x_min=None, x_max=None):
        if not self.obs_data:  # Ensure the file path is set for observed data
            self.status_msg = "ERROR: No file selected for Observed Data."
            self.status.setText(self.status_msg)
            self.status.setStyleSheet("color: red; font-weight: bold;") #set error message to red and bold
            return
        # if not self.diffx_df:  # Ensure the file path is set for observed data
        #     self.status_msg = "ERROR: No simulation ran to create Diffraction graph."
        #     self.status.setText(self.status_msg)
        #     self.status.setStyleSheet("color: red; font-weight: bold;") #set error message to red and bold
        #     return

        # Load the data from the .xy file (assuming the file is space or tab-separated)
        data = np.loadtxt(self.obs_data)
        
        # Extract x and y values (assuming the file has two columns)
        self.obs_df = self.obs_df.sort_values(by=["x"], ascending=True) #make sure x is non-decreasing order
        x = self.obs_df["x"]
        y = self.obs_df["y"]
        # Extracct x and y values for diffraction
        # TODO: update variable names to more appropriate ones later
        x2 = self.diffx_df["x"]
        y2 = self.diffx_df["y"]

         # Apply optional x-axis filtering
        if x_min is not None:
            mask = x >= x_min
            x, y = x[mask], y[mask]
            mask2 = x2 >= x_min
            x2, y2 = x2[mask2], y2[mask2]
        
        if x_max is not None:
            mask = x <= x_max
            x, y = x[mask], y[mask]
            mask2 = x2 <= x_max
            x2, y2 = x2[mask2], y2[mask2]
        
        # Create a Matplotlib figure and axis
        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot the x and y values as a line plot
        # Used selected colors for the lines
        self.obs_color_selected = self.obs_color.currentText().lower() 
        self.diffx_color_selected = self.diffx_color.currentText().lower()
        
        ax.plot(x, y, label='XY Data', color=self.obs_color_selected) #observe data
        ax.scatter(x2, y2, label='Diffraction Simulation Data', color=self.diffx_color_selected, s=5) #diffraction data
        
        # Label the graph
        ax.set_title("Graph of XY Data")
        ax.set_xlabel("Q (A^-1)") #This is currently 2*theta -> TODO: need to change to Q (should be around 1-10)
        # TODO: _init.py file hass tt_to_q() method to change from 2Theta to Q
        ax.set_ylabel("Intensity")
        # TODO: normalize intensity

        # Create a canvas and toolbar to embed the plot in the PyQt GUI
        self.canvas = FigureCanvas(fig)
        self.toolbar = NavigationToolbar(self.canvas, self)  
        
        # Add the canvas to the grid layout at the specified position (2, 9, 7, 11)
        self.grid.addWidget(self.toolbar, 14, 10, 1, 6)  
        self.grid.addWidget(self.canvas, 2, 9, 7, 11)
        self.canvas.draw()  # Render the canvas (display the plot)
        
        # set to true to indicate that the graph has been made
        self.is_graphed = True
 
    def update_graphs(self):
        # check if the files have been generated first
        if not self.is_graphed: 
            self.status_msg = "ERROR: Generate files first before updating graphs."
            self.status.setText(self.status_msg)
            self.status.setStyleSheet("color: red; font-weight: bold;")
        # check that x min and x max are not empty
        # elif not self.x_min.text() or not self.x_max.text():
        #     self.status_msg = "ERROR: Please input values for X Min and X Max to update graphs."
        #     self.status.setText(self.status_msg)
        #     self.status.setStyleSheet("color: red; font-weight: bold;")
        else:
            # update the graph with new values
            # get the new values from the text boxes
            # check that x_min and x_max are not empty
            if (self.x_min.text() and self.x_max.text()):
                x_min = float(self.x_min.text())
                x_max = float(self.x_max.text())
                # update the graph with new values
                self.create_graph(x_min, x_max)
            elif (self.x_min.text()):
                x_min = float(self.x_min.text())
                self.create_graph(x_min)
            elif (self.x_max.text()):
                x_max = float(self.x_max.text())
                self.create_graph(None, x_max)
            else:
                self.create_graph()
            # update status message
            self.status_msg = "Graphs Updated"
            self.status.setText(self.status_msg)
            self.status.setStyleSheet("color: green; font-weight: bold;")
            
        


def start_application():
    # Q application instance
    app = QApplication(sys.argv) #passes cmd line arguments
    # create window 
    window = MainWindow()
    window.showMaximized() #shows window
    # start event loop
    app.exec_()


#MAIN METHOD
start_application()


# -- FUTURE IMPLEMENTATION --
# include some sort of scrolling function to zoom into the graphs
# include a Transform button in order to adjust the numbers to different positions depending on change to x,y,z and alpha,gamma,beta
      
        
 
        # left layout
        # left_layout = QVBoxLayout()

        # # SECTION 1: INITIAL INFORMATION

        # # Create a horizontal layout for three input boxes
        # input_layout = QHBoxLayout()
        # # Create a label for "Wavelength"
        # wavelength_label = QLabel("Wavelength:")
        # input_layout.addWidget(wavelength_label)
        # for i in range(3):
        #     input_box = QLineEdit(self)  # Create input box
        #     input_box.setPlaceholderText(f"Input {i+1}")  # Optional placeholder text
        #     input_layout.addWidget(input_box)
        #     self.wavelength_inputs.append(input_box)  # Store the input box reference
  
        # # Create left layout position
        # left_layout.addLayout(input_layout)
        
        # # GENERATE BUTTON
        # # Generate button to generate .flts file  
        # create_button = QPushButton(self)
        # create_button.setText("Generate File")
        # # self.setCentralWidget(create_button)
        # create_button.clicked.connect(self.generate)
        # left_layout.addWidget(create_button, alignment=Qt.AlignHCenter)
        
        # # left panel widget
        # left_panel = QWidget()
        # left_panel.setLayout(left_layout)
        
        # # right panel widget
        # # empty for the time being
        # right_panel = QWidget()


        # # Splitter to hold left and right panels
        # splitter = QSplitter(Qt.Horizontal)
        # splitter.addWidget(left_panel)
        # splitter.addWidget(right_panel)
        # self.setCentralWidget(splitter)
        # splitter.setHandleWidth(3) # create vertical line
        # splitter.setStyleSheet("QSplitter::handle { background-color: gray; }")
        # splitter.setStretchFactor(0, 1)
        # splitter.setStretchFactor(1, 2)
        