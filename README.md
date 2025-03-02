# Human-Bioactives-Structure-Database

## 1.	INTRODUCTION

In silico drug discovery refers to the use of computational methods and technologies to identify, design and optimize potential drug candidates [1]. This approach involves the application of various techniques such as high throughput screening, virtual screening, molecular modelling, and so on to expedite the drug discovery process [2]. It offers several key benefits, such as predict toxicity, assess pharmacokinetics, and expedite the identification of promising compounds. 

These discovery processes rely heavily on pre-existing chemical database which provide a vast repository of molecular structures for virtual screening, lead optimization and drug design. The production of these databases is essential due to the exponential growth of biological and chemical knowledge, which has made it impractical to manually curate and organize the vast amount of data. These databases integrate information from various sources and provide a unified, searchable resource that links both chemical structure and biological relevance of the potential drug bio actives. They also help in the visualization of compounds to facilitate comparison of drug similarity and in silico drug target discovery [3]. 

## 2.	METHODS

### 2.1	DATABASE GENERATION AND ARCHITECHTURE

The Human Bioactives Structure Database is a chemi-informatics based database which provides information about the structure and molecular function of the bio actives. The database was generated to view 100 bio actives provided to us in a .sdf format. RDkit [4], which is an open source chemi-informatics library in python to get all the useful information from the raw file provided to us. A SQL database was build using SQLite [5] which was populated from the information obtained from the input file. Further, Tkinter [6], a library in python built for Interface production and visualization. Figure 1 shows an overview flowchart of the database production process.

<img width="1000" alt="Screenshot 2025-03-02 at 8 12 01 pm" src="https://github.com/user-attachments/assets/a9eaffba-48e2-41b4-ac52-5dc1bca85ac5" />

Figure 1: Flowchart of the database production process

### 2.2	DATABASE GENERATION

A set of Database Definition Language (DDL) statements were written through SQLite to create a database (see Protein_Databases.sql). Sqlite3 [7] library in python was used to create the database and insert data into the database. The DDL statements were loaded in python using sqlite3 and a database was created.A raw file in .sdf format was also given as input and iterated upon using RDkit. A Structured Data File (SDF) format file is a chemical data file format which contains information regarding the chemical structure and the associated data of the compounds. The raw file was iterated through, and the following information were extracted:

1.	Molecular name: “mol.GetProp(_name)” function was used to extract the name of the compound
2.	Chemical Structure: “AllChem.Compute2DCoords(mol)” function was used to obtain the coordinates of the atoms mentioned in the raw file. Then “Draw.MolToImage(mol)” function was used to convert the coordinates into a 2-D structure. The structure image was saved as a .png format in the output path. The image is then imported back into python and iterated over using libraries base64 [9] and PIL.Image [10] to convert the .png image into base64 encode by using “base64.b64encode(img_data).decode('utf-8')” function.
3.	SMILES: SMILES stands for Simplified Molecular Input Line Entry System which is used to translate a 3-D compound structure into a string of symbols which can be easily understood by various computer software [11]. “Chem.MolToSmiles(mol)” function was used to convert the coordinated present in the sdf file into a SMILES string.
4.	Chemical Fomula: “mol.GetProp('Formula')” function was used to extract the chemical formula of the compounds from tbe raw file.
5.	Molecular weight: “Descriptors.MolWt(mol)” function was used to calculate the molecular weight of the compunds.
6.	LogP: “Descriptors.MolLogP(mol)” function was used to calculate the LogP value of the compounds in the rawfile.
7.	LogD:  “mol.GetProp('LogD')” function was used to extract the LogD values from the raw file
8.	Number of atoms: “mol.GetNumAtoms()” function was used to obtain the number of atoms in the rawfile.
9.	Rotatable bonds: “Descriptors.NumRotatableBonds(mol)” was used to calculate the number of rotatable bonds for each compound in raw file.
10.	Ring Count: “Descriptors.RingCount(mol)” function was used to calculate the number of rings in a compound for each compound in the rawfile.
11.	Donors: “Descriptors.NumHDonors(mol)” function was used to calculate the number of H+ donors present in each compound in the raw file.
12.	Acceptors: “Descriptors.NumHAcceptors(mol)” function was used to calculate the number of H+ acceptors for each compound in the rawfile.
13.	Valence Electrons: “Descriptors.NumValenceElectrons(mol)” function was used to calculate the number of valences electrons present on each compound in the raw file.
14.	Polar Surface Area: Total Polar Surface Area (TPSA), which is defined as the surface sum over all polar atoms or molecules, primarily oxygen and nitrogen, including their attached hydrogen atoms is a widely used criteria for chemists to optimize a drug’s ability to permeate cells [12]. “Descriptors.TPSA(mol)” was used to calculate the TPSA of each compound in raw file.
15.	Fused Aromatic Rings: Fused aromatic rings (FAR) are a common structural motif in organic compounds, and their count can impact various properties such as aqueous solubility, lipophilicity, and the potential for drug-likeness [13]. Rings that are fused with each other share electrons as well. Therefore, to calculate FARs, we check for at least 2 electrons that are being shared between two rings. To do this, “Chem.GetSymmSSSR(mol)” was used to get the smallest set a symmetrical ring. Then a nested loop was run to get the length of both the rings. Futher, aromatic bonds were validated for each rings by using “mol.GetAtomWithIdx(atom).GetIsAromatic()”  function. If the rings are aromatic, they are check for any shared electrons between them. Every time, the number of shared atoms is greater than 1, the FAR count is increased to the number of rings fused together.

All these parameters were inferred from the raw file and inserted into the SQLite database using a Database Manipulation Language (DML) statement which used the keyword INSERT to command SQLite to insert the data into the respective columns in the database. Alongwith it a logger has been set to log and function as a check point for the success of each step as well as any errors that occurred during each step (6. Appendix , Figure 1). 



### 2.3	GRAPHICAL USER INTERFACE (GUI) PRODUCTION

TKinter was used to develop an interface for the database created in section 2.2. The entire GUI was divided into 3 frames to facilitate the layout of all the elements of the GUI. An overview and the frame placements of the three frames have been shown in Figure 3.


 <img width="1000" alt="image" src="https://github.com/user-attachments/assets/f35975bb-b7eb-44dc-9c98-043107564e7c" />


Figure 3: Overview of the GUI. The 3 frames have been shown as Frame 1, Fram 2 and Frame 3 as above.


The frames have been divided in the following manner:

1.	Frame 1: This frame incorporates the following:
a.	Human Bioactive Structure Database label created using the “Label()” in TKinter library.
b.	A search bar to search using “Entry()” function which allows user to type in a text space This input is stored in a variable which can be used for developing a function to search the data according to the molecular names.
c.	A submit button which is connected to a function “SearchDatabase()” that runs a SQL DML statement to search data database of the user input molecular name.
d.	A reset button which is connected to a “resetData()” function that clears the tree and loads the entire database back into the GUI. 
e.	A textbox showing the number of rows in currently in the database.
f.	Two drop-menus to sort the data into ascending and descending order. The drop-menus were created using the “Optionmenu()” function. The user can select the column which they wish to sort and the order in which they wish to sort.
g.	A sort button which is connect to the “sort()” function which takes in the values selected in the two menus and sort them accordingly. 
h.	An “export to .csv” button to allow the user to download data through the GUI. This button is connected to an “export()” function on the backend which writes the data into .csv format using csv library [14].
2.	Frame 2: This frame contains a treeview to display the main table. The function “displayData()” is called upon which call the database and fetched all the rows in the table using an SQL statement, converts base64 image into a visual 2-D structure image and inserts everything into treeview.
3.	Frame 3: This frame consists of all the filters and customizable parameter options mentioned as follows:
a.	A button to filter according to the Lipinski Rule of 5 criteria. The button is connected to a “FilterLipinskiRuleOfFive()” function that runs an SQL statement to filter the database according to the Lipinski Rule of 5 criteria.
b.	A button to filter according to lead likeness [16]. The button is connected to a “FilterLeadLikeness()” function which runs an SQL statement to filter the database according to the criteria for Lead likeness. 
c.	A button to filter according to the bioavailability criteria. The button is connected to a function “FilterBioAvailability()” that iterates through each row and checks how many criterias are satiated. It then adds +1 to the score for each passed criterion. If the score is greater an or equals to 6, the compound is displayed onto treeview.
d.	A clear filter button connected to “resetData()” to reset the treeview database.
e.	A number of customisable parameters as shown in Figure 3 which can be customised by the user. These parameters are connected to a customize button which runs a function “apply_paramters()” and dynamically builds an SQL query. The function checks whether both the textboxes are empty or have any value. If there is a value in any one of the check boxes, it further checks which text bos is empty for and switches the other with the default value already set according to the database. This is the appended into a dictionary. Then the length of the dictionary is measured and an sql statement will be built dynamically.
f.	A clear button connect to a “resetData()” function which resets the data back to original.

Apart from this, a logger in backend and message boxes on the frontend have been set to log and reports to the user, success of each step as well as any errors that occurred during each step. An image of the log as well as the message box have been mentioned in 6. Appendix, Figure 1 and Figure, respectively.


## 3.	RESULTS

### 3.1 HOW TO BROWSE THROUGH THE DATABASE

The homepage (Figure 3) is a common entry point for the GUI and showcases all the filters and parameter options available to use for the database along with the database and all the columns. The user can scroll through the treeview and have a look at all the columns (mentioned in section 2.2). The number of rows can be seen besides the rows label. This number changes according to the rows that have been filtered by a query. The search bar allows user to search according to molecular name by typing either partially or completely. The reset button resets the query and shows the entire database (Figure 3).


 <img width="1000" alt="image" src="https://github.com/user-attachments/assets/621678b7-04d1-43c3-af3f-b7b5c112ce73" />


Figure 4: Human Bioactive Structure DB: (A) Partial search by typing “dieth”. (B) Complete search by typing “dienestrol”.


The database can also be sorted according to the columns by selecting the column name from the drop-down menu and the order in which to sort from the ascending or descending drop-down menu (Figure 5). 

<img width="1000" alt="image" src="https://github.com/user-attachments/assets/0bbed8f2-cd4a-41b5-9328-c1a67069d66c" />

Figure 5: Human Bioactive Structure Database sorted according to Molecular weight in ascending order.

### 3.2 HOW TO QUERY THE DATABASE

There are three filter buttons under the “Filter Options” label allows you to filter the data based on the pre-specified criteria. The criteria shows up in the respective parameter boxes to make the user aware of the criteria (Figure 6).

 
<img width="1000" alt="image" src="https://github.com/user-attachments/assets/84ac981d-7185-4f29-9b99-7efc0ecbecaf" />

Figure 6: Human Bioactive Structure Database filtered by Lipinski Rule of 5 button.

These filters can also be customized by adding or deleting the as many criterions according to the user (Figure 7). 

 <img width="1000" alt="image" src="https://github.com/user-attachments/assets/e3b66469-fa7c-4a41-bc1e-4bbc377d21cb" />

Figure 7: Human Bioactive Structure Database filtered by user specifies parameters

### 3.3 HOW TO DOWNLOAD DATA 

The “export to .csv” button allows the user to save the queried database as a .csv format. Upon pressing the button a dialog box appears where the user can enter the file name and the path where they want to save the csv file (Figure 8). The file format has been mentioned in 6. Appendix, Figure 3.

<img width="1000" alt="image" src="https://github.com/user-attachments/assets/11b8976b-e87c-4c77-9686-a4ac734b2386" />

Figure 8: Human Bioactive Structure Database dialog box asking for user specifies filename and path to save the data


## 4.	IMPROVEMENTS


The following improvements can be implemented in the database:

i.	An autofill search bar that allows to facilitate searching of data.
ii.	To enable sort function to sort in any order for the different filter options and customizable parameters.


## 5.	REFERENCES

1.	Shaker, B., Ahmad, S., Lee, J., Jung, C. and Na, D. (2021). In silico methods and tools for drug discovery. Computers in Biology and Medicine, 137, p.104851. doi:https://doi.org/10.1016/j.compbiomed.2021.104851.
2.	Ekins, S., Mestres, J. and Testa, B. (2007). In silicopharmacology for drug discovery: methods for virtual ligand screening and profiling. British Journal of Pharmacology, [online] 152(1), pp.9–20. doi:https://doi.org/10.1038/sj.bjp.0707305.
3.	Kim, C. and Kim, E. (2019). Rational Drug Design Approach of Receptor Tyrosine Kinase Type III Inhibitors. Current Medicinal Chemistry, [online] 26(42), pp.7623–7640. doi:https://doi.org/10.2174/0929867325666180622143548.
4.	www.rdkit.org. (n.d.). RDKit. [online] Available at: https://www.rdkit.org.
5.	SQLite (2019). SQLite Home Page. [online] Sqlite.org. Available at: https://www.sqlite.org/index.html.
6.	Python Software Foundation (2019). tkinter — Python interface to Tcl/Tk — Python 3.7.2 documentation. [online] python.org. Available at: https://docs.python.org/3/library/tkinter.html.
7.	Python Software Foundation (n.d.). sqlite3 — DB-API 2.0 interface for SQLite databases — Python 3.8.2 documentation. [online] docs.python.org. Available at: https://docs.python.org/3/library/sqlite3.html.
8.	lifechemicals.com. (n.d.). How to work with Structured Data Files (SDF files) | Order and Supply | Life Chemicals. [online] Available at: https://lifechemicals.com/order-and-supply/how-to-work-with-sd-files [Accessed 18 Feb. 2024].
9.	docs.python.org. (n.d.). base64 — Base16, Base32, Base64, Base85 Data Encodings — Python 3.9.5 documentation. [online] Available at: https://docs.python.org/3/library/base64.html.
10.	pillow.readthedocs.io. (n.d.). Pillow. [online] Available at: https://pillow.readthedocs.io/en/stable/index.html.
11.	Sustainable Futures / P2 Framework Manual 2012 EPA-748-B12-001 Appendix F. SMILES Notation Tutorial. (n.d.). Available at: https://www.epa.gov/sites/default/files/2015-05/documents/appendf.pdf.
12.	Barret, R. (2018). Importance and Evaluation of the Polar Surface Area (PSA and TPSA). Therapeutical Chemistry, pp.89–95. doi:https://doi.org/10.1016/b978-1-78548-288-5.50005-6.
13.	Ritchie, T.J. and Macdonald, S.J.F. (2009). The impact of aromatic ring count on compound developability – are too many aromatic rings a liability in drug design? Drug Discovery Today, 14(21-22), pp.1011–1020. doi:https://doi.org/10.1016/j.drudis.2009.07.014.
14.	Python Software Foundation (2020). csv — CSV File Reading and Writing — Python 3.8.1 documentation. [online] Python.org. Available at: https://docs.python.org/3/library/csv.html.
15.	Lipinski, C.A., Lombardo, F., Dominy, B.W. and Feeney, P.J. (1997). Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. Advanced Drug Delivery Reviews, [online] 23(1-3), pp.3–25. doi:https://doi.org/10.1016/s0169-409x(96)00423-1.
16.	Congreve, M., Carr, R., Murray, C. and Jhoti, H. (2003). A ‘Rule of Three’ for fragment-based lead discovery?. Drug Discovery Today, 8(19), pp.876–877. doi:https://doi.org/10.1016/s1359-6446(03)02831-9.
