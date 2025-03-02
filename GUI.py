#Import all modules
#from 2743102_Database import * ##  Uncomment this to run both the scripts at once.
import logging
from tkinter import *
from tkinter import ttk
from tkinter import messagebox,filedialog
import pandas as pd
import csv
import os
import sqlite3
from PIL import Image, ImageTk
import io
import base64
import numpy as np 
import sys
from io import BytesIO
from datetime import datetime

#Setting up the logger###
logger = logging.getLogger()
logger.setLevel(logging.INFO)
fh_log = logging.FileHandler('2743102_log.txt')
fh_log.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(fh_log)

#Background colour###
bkg = '#da8d9b'
bkg1 = '#e9ecee'
db_file = '/Users/harshitasrivastava/Library/Mobile Documents/com~apple~CloudDocs/UofG Bioinfo/Semester 2/Chemical Structures/Chemical Structures /Chemical Databases/protein.db'

###INITIALIZE GLOBAL PARAMETERS###
#Global list to store image references
image_list = []

#Global list to store query rows 
query_rows = []

#Setting defaults for the customisable parameters
molwt_min, molwt_max = '124.23','807.89'
logp_min,logp_max = '-7.95','12.03'
logd_min, logd_max = '-17.9','9.64'
rc_min, rc_max = '0','10'
rb_min,rb_max = '0','39'
d_min,d_max = '0','15'
a_min,a_max = '0','15'
tpsa_min,tpsa_max = '0','334.59'
far_min,far_max ='0','6'

###Connect to Database###
conn = sqlite3.connect(db_file)
#Create a cursor
cur = conn.cursor()

###configure workspace###
ws = Tk()
ws.title("Human Bioactive StructureDB")
ws.geometry('1366x1080')
ws.configure(bg=bkg)

###FUNCTIONS###
#Create a function to insert values into treeview
def insertDataToTreeview(rows):
    try:
        global query_rows
        query_rows = []
        #Loop through the list to insert the data into the GUI.
        count=0
        #Iterate over each row and add it to tree view
        for row in rows:
            bytes_decoded = base64.b64decode(row[1])
            img = Image.open(BytesIO(bytes_decoded))
            out_jpg = img.convert("RGB")
            out_jpg = out_jpg.resize((100, 100))
            img_import = ImageTk.PhotoImage(out_jpg)
            image_list.append(img_import)
            trv.insert(parent='', index='end',image= img_import,values=(row[0],row[2],row[3],row[4],row[5],row[6],row[7],row[8],row[9],row[10],row[11],row[12],row[13],row[14]))
            #Set the fetched data into variable globally accessible by other functions
            query_rows.append([row[0],row[2],row[3],row[4],row[5],row[6],row[7],row[8],row[9],row[10],row[11],row[12],row[13],row[14]])
            count+=1
        row_count.set(count)
    except Exception as e:
        logger.error(f'Exception occured while inserting rows from Database into GUI: {e}')
        messagebox.showerror('ERROR',f"An error occurred while retrieving data:\n{str(e)}")
        pass
    else:
        logger.info(f'Rows inserted!')


#Create a function to export data to a CSV file
def export():
    try:
        global query_rows
        #If there are  no entries in the table, show an error otherwise open a dialog box to save the csv file
        if len(query_rows) < 1:
            messagebox.showerror("No Data", " No data available to export")
            return False
        fln = filedialog.asksaveasfilename(initialdir=os.getcwd(), title="Save CSV", filetypes=(("CSV File", "*.csv"),("All Files","*.*")))
        logger.info("Export path = {}".format(fln))

        #Write the output into the file
        with open(fln,'w') as myfile:
            exp_writer = csv.writer(myfile, delimiter = ',')
            myfile.write(f'Molecular_Name, SMILE,Chemical_Formula,Molecular_Weight,LogP,LogD,Number_of_atoms,Rotatable_Bonds,Number_of_rings,Donors,Acceptors,Valence_Electrons,Polar_Surface_Area,Fused_Aromatic_Rings\n')
            for i in query_rows:
                exp_writer.writerow(i)
    except Exception as e:
        logger.error(f'An exception occured while exporting file as CSV: {e}')
        pass
    else:
        logger.info("File saved successfully. Path = {}".format(os.path.basename(fln)))
        messagebox.showinfo(title="Data Exported Successfully!",message="Data Exported Successfully!. Path = {}".format(os.path.basename(fln)))

#Function to reset the database
def resetData():
    try:
        #Reinitialize Customizable parameters
        molwt_text1.set(''),molwt_text2.set(''),logp_text1.set(''),logp_text2.set(''),logd_text1.set(''),logd_text2.set(''),rc_text1.set(''),rc_text2.set(''),d_text1.set(''),d_text2.set(''),a_text1.set(''),a_text2.set(''),rb_text1.set(''),rb_text2.set(''),tpsa_text1.set(''),tpsa_text2.set(''),far_text1.set(''),far_text2.set('')

        #Clear the Treeview
        for record in trv.get_children():
            trv.delete(record)

        #Fetch all the records from SQL database
        cur.execute("SELECT * FROM PROTEIN_DATABASE;")
        rows = cur.fetchall()

        #Loop through the list to insert the data into the GUI.
        #Clear the image list     
        image_list.clear()
        #Set a counter to count the number of rows
        count =0
        #Iterate over each row and add it to tree view
        for row in rows:
            bytes_decoded = base64.b64decode(row[1])
            img = Image.open(BytesIO(bytes_decoded))
            out_jpg = img.convert("RGB")
            out_jpg = out_jpg.resize((100, 100))
            img_import = ImageTk.PhotoImage(out_jpg)
            image_list.append(img_import)
            trv.insert(parent='', index='end',image= img_import,values=(row[0],row[2],row[3],row[4],row[5],row[6],row[7],row[8],row[9],row[10],row[11],row[12],row[13],row[14]))
            count+=1
        row_count.set(count)
    except Exception as e:
        logger.error(f"Error occured while clearing the Database:{e}")
        messagebox.showerror('Error','Failed to Clear the Database')
    else:
        logger.info(f'The Database has been Cleared Successfully!')
    

#Create a function to display the data into treeview
def displayData():
    try:
        #Fetch all the records from SQL database
        cur.execute("SELECT * FROM PROTEIN_DATABASE;")
        rows = cur.fetchall()
        #Call the function to insert the data into TreeView
        insertDataToTreeview(rows)
    except Exception as e:
        logger.error(f'An exception while displaying data on TreeView: {e}')
        pass
    else:
        logger.info(f"Displaying data from Database Successful.")

#Create a function to apply customisable parameters
def apply_parameters():
    try:
        #Clear the Treeview
        for record in trv.get_children():
            trv.delete(record)

        #Initialize empty dictionary
        filter_dict = {}
        
        #Check which checkboxes are selected and set their values in the dictionary accordingly
        #Molecular Weight
        if molwt_text1.get()+molwt_text2.get() != '':
            if molwt_text1.get() == '':
                filter_dict['molwt'] = [molwt_min,molwt_text2.get()]
            elif molwt_text2.get() == '':
                filter_dict['molwt']=[molwt_text1.get(),molwt_max]
            else:
                filter_dict['molwt']=[molwt_text1.get(), molwt_text2.get()]
    
        #LogP
        if logp_text1.get()+logp_text2.get() != '':
            if logp_text1.get() == '':
                filter_dict['logP'] = [logp_min,logp_text2.get()]
            elif logp_text2.get() == '':
                filter_dict['logP']=[logp_text1.get(),logp_max]
            else:
                filter_dict['logP']=[logp_text1.get(), logp_text2.get()]

        #LogD
        if logd_text1.get()+logd_text2.get() != '':
            if logd_text1.get() == '':
                filter_dict['logD'] = [logd_min,logd_text2.get()]
            elif logd_text2.get() == '':
                filter_dict['logD']=[logd_text1.get(),logd_max]
            else:
                filter_dict['logD']=[logd_text1.get(), logd_text2.get()]

        #Ring Count
        if rc_text1.get()+rc_text2.get() != '':
            if rc_text1.get() == '':
                filter_dict['ringcount'] = [rc_min,rc_text2.get()]
            elif rc_text2.get() == '':
                filter_dict['ringcount']=[rc_text1.get(),rc_max]
            else:
                filter_dict['ringcount']=[rc_text1.get(), rc_text2.get()]

        #Donors
        if d_text1.get()+d_text2.get() != '':
            if d_text1.get() == '':
                filter_dict['H_donors'] = [d_min,d_text2.get()]
            elif d_text2.get() == '':
                filter_dict['H_donors']=[d_text1.get(),d_max]
            else:
                filter_dict['H_donors']=[d_text1.get(), d_text2.get()]    
        
        #Acceptors
        if a_text1.get()+a_text2.get() != '':
            if a_text1.get() == '':
                filter_dict['H_acceptors'] = [a_min,a_text2.get()]
            elif logp_text2.get() == '':
                filter_dict['H_acceptors']=[a_text1.get(),a_max]
            else:
                filter_dict['H_acceptors']=[a_text1.get(),a_text2.get()]
        
        #Rotatable Bonds
        if rb_text1.get()+rb_text2.get() != '':
            if rb_text1.get() == '':
                filter_dict['rotatablebonds'] = [rb_min,rb_text2.get()]
            elif rb_text2.get() == '':
                filter_dict['rotatablebonds']=[rb_text1.get(),rb_max]
            else:
                filter_dict['rotatablebonds']=[rb_text1.get(),rb_text2.get()]
        
        #Polar Surface Area(PSA)
        if tpsa_text1.get()+tpsa_text2.get() != '':
            if tpsa_text1.get() == '':
                filter_dict['polarSA'] = [tpsa_min,tpsa_text2.get()]
            elif tpsa_text2.get() == '':
                filter_dict['polarSA']=[tpsa_text1.get(),tpsa_max]
            else:
                filter_dict['polarSA']=[tpsa_text1.get(), tpsa_text2.get()]
        
        #Fused Aromatic Ring
        if far_text1.get()+ far_text2.get() != '':
            if far_text1.get() == '':
                filter_dict[''] = [far_min,far_text2.get()]
            elif far_text2.get() == '':
                filter_dict['logP']=[far_text1.get(),far_max]
            else:
                filter_dict['logP']=[far_text1.get(), far_text2.get()]

        #Buiild the sql statement dynamically
        #Initialise
        sql_statement = ''
        #Check the length of the sql statement and execute accordingly
        if len(filter_dict)==1:
            for key,value in filter_dict.items():
                sql_statement =  f'SELECT * FROM PROTEIN_DATABASE WHERE {key} BETWEEN {value[0]} AND {value[1]} ;'
                cur.execute(sql_statement)
                logger.info(f'Statement for customizable parameters: {sql_statement}')
        elif len(filter_dict)>1:
            sql_statement = 'SELECT * FROM PROTEIN_DATABASE WHERE '
            for key,value in filter_dict.items():
                sql_statement += f'{key} BETWEEN {value[0]} AND {value[1]} AND ' 
            sql = sql_statement[:-4]+ ';' #Remove the last AND and ; instead
            cur.execute(sql)
            logger.info(f'Statement for customizable parameters: {sql}')
        #Fetch all the rows from the Database and add them to the TreeView
        rows = cur.fetchall()
        #Call the function to insert the data into TreeView
        insertDataToTreeview(rows)
    except Exception as e:
        #logger.error(f"Error occurred while filtering the database:{e}")
        messagebox.showerror(title= "Error",message="An error occured while trying to apply filters")
        resetData()
        pass
    else:
        logger.info(f"Database filtered successfully acoording to user mentioned parameters")
        messagebox.showinfo(title= "Filtering Successful",message="The database has been filtered according to your specified criteria.")

#Add function to search database
def SearchDatabase():
    try:
        #Get the value of the searchbox
        lookup_record = entry.get()
        #Clear the Treeview
        for record in trv.get_children():
            trv.delete(record)
        sql = f"SELECT * FROM PROTEIN_DATABASE WHERE molname LIKE '%{lookup_record}%';"
        cur.execute(sql)
        rows = cur.fetchall()
        #Call the function to insert the data into TreeView
        insertDataToTreeview(rows)
    except Exception as e:
        logger.error(f"An error occurred while searching the Database:{e}")
        messagebox.showerror('Error', 'Search Failed! Please check your input and try again.')
    else:
        if len(rows) == 0:
            logger.info(f'No Results Found!', f'No records found with "{lookup_record}"')
            messagebox.showinfo('No Results','No records found!')
            resetData()
            pass
        else:
            logger.info(f'Found {len(rows)} Records for {lookup_record}!')
            messagebox.showinfo('Found Records','Records Found!')
        
#Create a function to sort the data
def sort():
    try:
        #Get the opeions selected by user
        sort_var = clicked.get()
        sort_order = clicked_order.get()

        #Initialize an empty string
        col_name = ''

        #Clear the Treeview
        for record in trv.get_children():
            trv.delete(record)

        #Fetch the rows by ORDER from the DATABASE
        sort_dict = {'Molecular Name': 'molname', 'Molecular Weight': 'molwt', 'LogP':'logP','Number of Rings':'ringcount','Donors':'H_donors','Acceptors':'H_acceptors','Polar Surface Area':'polarSA','Fused Aromatic Rings':'FAR'}

        #Iterate through dictionary and check if the key is selected, then add it to col_name with appropriate order
        for key in sort_dict.keys():
            if sort_var == key and sort_order == 'Ascending':
                col_name = sort_dict[key]
                sql_ascend =  f'SELECT * FROM PROTEIN_DATABASE ORDER BY {col_name} COLLATE NOCASE ASC;'
                cur.execute(sql_ascend)
            elif sort_var == key and sort_order == 'Descending':
                col_name = sort_dict[key]
                sql_descend = f"SELECT * FROM PROTEIN_DATABASE ORDER BY {col_name} COLLATE NOCASE DESC;"
                cur.execute(sql_descend)
        #Fetch all the rows from the Database and add them to the TreeView
        rows = cur.fetchall()
        #Call the function to insert the data into TreeView
        insertDataToTreeview(rows)
    except Exception as e:
        logger.error(f"An error occured while sorting the DataBase:{e}")
        messagebox.showerror('Error','Sorting Failed! Please Check Your Inputs and Try Again.')
        resetData()
        pass
    else:
        logger.info(f'The table has been sorted by {sort_var} in {sort_order} order.')
        messagebox.showinfo(message='SORT SUCCESSFUL!')
    

#Create a function to filter Lipinski rule of 5
def FilterLipinskiRuleOfFive():
    try:
        #Clear the Treeview
        for record in trv.get_children():
            trv.delete(record)
        sql = f'SELECT * FROM PROTEIN_DATABASE  WHERE  molwt<=500.00 AND logP <=5.0  AND H_donors <=5 AND H_acceptors <= 10;'
        cur.execute(sql)
        rows = cur.fetchall()
        #Call the function to insert the data into TreeView
        insertDataToTreeview(rows)
        #Set filter parameter values to customizble parameters
        molwt_text2.set(500)
        logp_text2.set(5.0)
        d_text2.set(5)
        a_text2.set(10)
    except Exception as e:
        logger.error(f"An error occurred while filtering using Lipinski Rule Of Five:{e}")
        messagebox.showerror(f'ERROR - Filtering Failed!{e}')
        resetData()
        pass
    else:
        logger.info(f'Filtered results based on Lipinski Rule of 5')
        messagebox.showinfo(title='Filtering Complete',message="Data has been filtered successfully.")

#Create a function to filter Lead Likeness
def  FilterLeadLikeness():
    try:
        #Clear the Treeview
        for record in trv.get_children():
            trv.delete(record)
        sql = f'SELECT * FROM PROTEIN_DATABASE WHERE molwt<=450 AND logD>=-4  AND logD<=4 AND ringcount<=4 AND rotatablebonds<=10   AND H_donors<=5 AND H_acceptors<=8;'
        cur.execute(sql)
        rows = cur.fetchall()
        #Set filter parameter values to customizble parameters
        molwt_text2.set(450)
        logd_text1.set(-4)
        logd_text2.set(4)
        rc_text2.set(4)
        rb_text2.set(10)
        #Call the function to insert the data into TreeView
        insertDataToTreeview(rows)
    except  Exception as e:
        logger.error(f"An error occurred while filtering using Lead Likeness:{e}")
        messagebox.showerror(f'ERROR - Filtering Failed! {e}')
        resetData()
    else:
        logger.info('Filtered results based on Lead Likeness')
        messagebox.showinfo(title='Filtering Complete',message="Data has been filtered successfully.")
        

#Create a function to calculate the bioavailability score
def bioavalability(db_file):
    try:
        #Call filter function to retrived all the rows that pass the Bioavailability criterion and check for their scores. 
        #If the score is >=6, the compound passes as bioavaliable
        sql = f'SELECT * FROM PROTEIN_DATABASE;'
        cur.execute(sql)
        all_rows = cur.fetchall()

        bioavailabilty_query_list = [] #To append the list of molecules that pass as bioavailable
        for i in range (0, len(all_rows)):
            score = 0
            list = all_rows[i]
            if list[4] <= 500.0: #Molecular wt
                score+=1
            if list[5] <= 5.0: #LogP
                score+=1
            if list[8] <=10: #Rotatble bonds
                score+=1
            if list[10] <= 5: #H_Donors
                score+=1
            if list[11]<=10: #H_Acceptors
                score+=1
            if list[13]<=200.0: #Polar Surface Area
                score+=1
            if list[14]<=5: #Fused Aromatic Ring
                score+=1
            #To check for bioavalability score
            if score >= 6: 
                logger.info(f'{list[0]} score = {score}')
                logger.info('Compound Passed The Filter Criterion.\n')
                bioavailabilty_query_list.append(all_rows[i])
        logger.info(f'Number of compunds that pass the bioavaliability criteria = {len(bioavailabilty_query_list)}.')
        logger.info(f'Bioavailability Score Calculation Completed Successfully!\n')
        return bioavailabilty_query_list
    except Exception as e:
        logger.error(f"Error occurred while calculating Bioavailability Score :\n{e}")
        pass
        

#Create a function to filter Bioavailability
def  FilterBioAvailability():
    try:
        #Clear the Treeview
        for record in trv.get_children():
            trv.delete(record)
        bioavailabilty_query_list = bioavalability(db_file)
        #Set filter parameter values to customizble parameters
        molwt_text2.set(500)
        logp_text2.set(5)
        d_text2.set(5)
        a_text2.set(10)
        rb_text2.set(10)
        tpsa_text2.set(200)
        far_text2.set(5)
        #Call the function to insert the data into TreeView
        insertDataToTreeview(bioavailabilty_query_list)
    except Exception as e:
        logger.error(f"Error Occurred While Filtering Data:\n{e}\n")
        messagebox.showerror(f'ERROR - Filtering Failed! {e}')
    else:
        logger.info(f'Data Filtered Successfully!\n')
        messagebox.showinfo(title='Filtering Complete',message="Data has been filtered successfully.")
        
##FRAME 1##
#Add a frame for the search bar
search_frame = Frame(ws, bg = bkg, height = 500, width = 1080, borderwidth=5)
search_frame.pack(side = 'top', expand= 'yes',padx=7, pady=7, anchor='center')

title = Label(search_frame, text="Human Bioactive Structure Database", font=('Times New Roman', 44,'bold'), fg='black', bg=bkg, anchor= 'e')
title.grid(row=0, column=1, padx=7, pady=7)

my_label = Label(search_frame, text="Search Bioactive names:", font=('Times New Roman', 20,'bold'), fg='black', bg=bkg)
my_label.grid(row=1, column=0, padx=7, pady=7)

#Create an entry box

entry = Entry(search_frame, font=('Times New Roman', 14), justify='left', bg='white', width = 140, fg='black')
entry.grid(row=1, column=1,  padx=7, pady=7)

#Add a search button
search_button = Button(search_frame, text = "Submit", font =('Times New Roman', 16,'bold'),bg= bkg ,command = SearchDatabase, width = 12)
search_button.grid(row=1, column=2, padx=7, pady=7)

#Add Reset Button
reset_button = Button(search_frame, text = "Reset", font = ('Times New Roman', 16,'bold'),bg=bkg,command = resetData, width=9)
reset_button.grid(row=1,column=3,padx=7,pady=7)

#Create a label for sorting the table
sort_label = Label(search_frame, text="Sort Data:", font=('Times New Roman', 19,'bold'), fg='black', bg=bkg)
sort_label.grid(row=4, column=1, sticky='E', padx= 7, pady= 7)   

#Create a dropdown menu for the sorting data
clicked = StringVar()
clicked.set("--SELECT--")
drop = OptionMenu(search_frame, clicked,"--SELECT--","Molecular Name","Molecular Weight","LogP","Number of Rings","Donors","Acceptors","Polar Surface Area","Fused Aromatic Rings")
drop.grid(row=4, column=2, sticky='E', padx= 7, pady= 7)

#Create a dropdown to select how the data will be selected
clicked_order = StringVar()
clicked_order.set("Ascending")
drop1 = OptionMenu(search_frame, clicked_order,"Ascending","Descending")
drop1.grid(row=4, column=3, sticky='E', padx= 7, pady= 7)

#Button to run the sorted value
sort_button = Button(search_frame, text = "Sort", font = ('Times New Roman', 16,'bold'),bg=bkg,command = sort)
sort_button.grid(row=4,column=4,padx=7,pady=7)

#Button to save the data as csv
export_button = Button(search_frame, text = "Export as .csv", font = ('Times New Roman', 16,'bold'),bg=bkg,command=export)
export_button.grid(row=4,column=0,padx=7,pady=7,sticky='W')

#Create a label for number of rows
row_label = Label(search_frame, text="Rows:", font=('Times New Roman', 19,'bold'), fg='black', bg=bkg)
row_label.grid(row=4,column=0,padx=7,pady=7,sticky='E')

#Create a textbox to display the number of rows
row_count = StringVar()
row_text = Entry(search_frame, font=('Times New Roman', 14),textvariable= row_count ,justify='left', width = 9, bg='white', fg='black')
row_text.grid(row=4,column=1,padx=7,pady=7,sticky='W')

##FRAME 2##
#Add a frame for Database
frame = Frame(ws, bg = bkg, width = 1300, height = 900, relief = GROOVE)
frame.pack_propagate(False) #Prevent the y
frame.pack(side ='right',pady = 10, padx =10)


#Create a style for the treeview
style = ttk.Style()
style.theme_use("clam")
style.configure(".", font=("Times New Roman", 18))
style.configure("Treeview.Heading", font=("Times New Roman", 16, "bold"))
style.configure("Treeview", rowheight=100) # set row height
style.configure("Treeview", width =100)    # set cell width

#Create Treeview Scrollbar
#Y axis scrollbar
scroll_y = ttk.Scrollbar(frame, orient = "vertical")
scroll_y.pack(side= RIGHT, fill= Y)

#X axis scrollbar
scroll_x = ttk.Scrollbar(frame, orient="horizontal")
scroll_x.pack(side=BOTTOM, fill=X)

#Create Treeview
trv = ttk.Treeview(frame, columns = (1,2,3,4,5,6,7,8,9,10,11,12,13,14),  yscrollcommand= scroll_y.set, xscrollcommand= scroll_x.set)

#Configure the scrollbar
scroll_y.config(command= trv.yview)
scroll_x.config(command =  trv.xview)

#Configure Columns
trv.column("#0",anchor=CENTER)
trv.column(1,anchor='w',width=300)
trv.column(2,anchor=CENTER,width=600)
trv.column(3,anchor=CENTER,width=150)
trv.column(4,anchor=CENTER,width=150)
trv.column(5,anchor=CENTER,width=80)
trv.column(6,anchor=CENTER,width=80)
trv.column(7,anchor=CENTER,width=150)
trv.column(8,anchor=CENTER,width=200)
trv.column(9,anchor=CENTER,width=150)
trv.column(10,anchor=CENTER,width=150)
trv.column(11,anchor=CENTER,width=150)
trv.column(12,anchor=CENTER,width=150)
trv.column(13,anchor=CENTER,width=150)
trv.column(14,anchor=CENTER,width=170)

#Add heading to columns
trv.heading("#0", text="Chemical Structure")
trv.heading(1, text="Molecular Name")
trv.heading(2, text="SMILES")
trv.heading(3,text="Chemical Formula")
trv.heading(4,text="Molecular Weight")
trv.heading(5,text="LogP")
trv.heading(6,text="LogD")
trv.heading(7,text="Number of atoms")
trv.heading(8,text="Rotatable Bond Number")
trv.heading(9,text="Number of rings")
trv.heading(10,text="Donors")
trv.heading(11,text="Acceptors")
trv.heading(12,text="Valence Electrons")
trv.heading(13,text="Polar Surface Area")
trv.heading(14,text="Fused Aromatic Rings")

#Call the displayData function to load the database into the GUI
displayData()

#Pack to screen
trv.pack(fill ='both',expand = 'yes')


##FRAME 3##
filter_frame = Frame(ws, bg = bkg1,width = 400, height= 851 ,borderwidth=5, relief = "groove")
filter_frame.pack_propagate(False)
filter_frame.pack(side = 'left', expand= 'yes', padx= 7)

#Create a label for filters
filter_label = Label(filter_frame, text="FILTER OPTIONS", font=('Times New Roman', 20,'bold'), fg='black', bg=bkg1)
filter_label.place(relx=.5, rely=.05, anchor='center')

#Add Lipinski filter Button
Lipinski_button = Button(filter_frame, text = "Lipinski Rule of 5", font =('Times New Roman', 20),bg= bkg1 ,width = 20 ,command = FilterLipinskiRuleOfFive)
Lipinski_button.place(relx=.5, rely=.1, anchor='center')

#Add Lead Likeness Filter Button
Leadlikeness_button = Button(filter_frame, text = "Lead Likeness", font =('Times New Roman', 20),bg= bkg1 ,width = 20 ,command = FilterLeadLikeness)
Leadlikeness_button.place(relx=.5, rely=.15, anchor='center')

#Add Bioavailability Filter Button
Bioavailability_button = Button(filter_frame, text = "Bioavailability", font =('Times New Roman', 20),bg= bkg1 ,width = 20 ,command = FilterBioAvailability)
Bioavailability_button.place(relx=.5, rely=.2, anchor='center')

#Add a Clear Filter button
Clearfilter_button = Button(filter_frame, text = "Clear Filter", font =('Times New Roman', 20),bg= bkg1 ,width = 20 ,command = resetData)
Clearfilter_button.place(relx=.5, rely=.3, anchor='center')

#Add label for Customizable Parameters
parameter_label = Label(filter_frame, text="PARAMETERS", font=('Times New Roman', 20,'bold'), fg='black', bg=bkg1)
parameter_label.place(relx=.5, rely=.40, anchor='center')

subparameter_label = Label(filter_frame, text="(Enter Numbers Only)", font=('Times New Roman', 15,'bold'), fg='black', bg=bkg1)
subparameter_label.place(relx=.5, rely=.43, anchor='center')

from_label = Label(filter_frame, text="FROM", font=('Times New Roman', 14,'bold'), fg='black', bg=bkg1)
from_label.place(relx=.52, rely=0.47, anchor='w')

to_label = Label(filter_frame, text="TO", font=('Times New Roman', 14,'bold'), fg='black', bg=bkg1)
to_label.place(relx=.73, rely=0.47, anchor='w')

#Add Customizable parameters
#Molecular Weight
molwt_label = Label(filter_frame, text="Molecular Weight", font=('Times New Roman', 17,'bold'), fg='black', bg=bkg1)
molwt_label.place(relx=.51, rely=0.5, anchor='e')

molwt_text1,molwt_text2 = StringVar(),StringVar()
mw_text1 = Entry(filter_frame, font=('Times New Roman', 14),textvariable= molwt_text1 ,justify='left', width = 9, bg='white', fg='black')
mw_text1.place(relx=.52, rely=0.5, anchor='w')

mw_text2 = Entry(filter_frame, font=('Times New Roman', 14),textvariable= molwt_text2 ,justify='left', width= 9,  bg='white', fg='black')
mw_text2.place(relx=.73, rely=0.5, anchor='w')

#LogP
lp_label = Label(filter_frame, text="LogP", font=('Times New Roman', 17,'bold'), fg='black', bg=bkg1)
lp_label.place(relx=.51, rely=0.55, anchor='e')

logp_text1,logp_text2 = StringVar(),StringVar()
lp_text1 = Entry(filter_frame, font=('Times New Roman', 14),textvariable= logp_text1 ,justify='left', width = 9, bg='white', fg='black')
lp_text1.place(relx=.52, rely=0.55, anchor='w')

lp_text2 = Entry(filter_frame, font=('Times New Roman', 14),textvariable= logp_text2 ,justify='left', width= 9,  bg='white', fg='black')
lp_text2.place(relx=.73, rely=0.55, anchor='w')

#LogD
ld_label = Label(filter_frame, text="LogD", font=('Times New Roman', 17,'bold'), fg='black', bg=bkg1)
ld_label.place(relx=.51, rely=0.60, anchor='e')

logd_text1,logd_text2 =StringVar(),StringVar()
ld_text1 = Entry(filter_frame, font=('Times New Roman', 14),textvariable=logd_text1 ,justify='left', width = 9, bg='white', fg='black')
ld_text1.place(relx=.52, rely=0.60, anchor='w')

ld_text2 = Entry(filter_frame, font=('Times New Roman', 14),textvariable= logd_text2 ,justify='left', width= 9,  bg='white', fg='black')
ld_text2.place(relx=.73, rely=0.60, anchor='w')

#Ring Count
rcount_label = Label(filter_frame, text="Ring Count", font=('Times New Roman', 17,'bold'), fg='black', bg=bkg1)
rcount_label.place(relx=.51, rely=0.65, anchor='e')

rc_text1,rc_text2 = StringVar(),StringVar()
rcount_text1 = Entry(filter_frame, font=('Times New Roman', 14),textvariable= rc_text1 ,justify='left', width = 9, bg='white', fg='black')
rcount_text1.place(relx=.52, rely=0.65, anchor='w')

rcount_text2 = Entry(filter_frame, font=('Times New Roman', 14),textvariable= rc_text2 ,justify='left', width= 9,  bg='white', fg='black')
rcount_text2.place(relx=.73, rely=0.65, anchor='w')

#Donors
donor_label = Label(filter_frame, text="Donors", font=('Times New Roman', 17,'bold'), fg='black', bg=bkg1)
donor_label.place(relx=.51, rely=0.70, anchor='e')

d_text1,d_text2 = StringVar(),StringVar()
donor_text1 = Entry(filter_frame, font=('Times New Roman', 14),textvariable=d_text1 ,justify='left', width = 9, bg='white', fg='black')
donor_text1.place(relx=.52, rely=0.70, anchor='w')

donor_text2 = Entry(filter_frame, font=('Times New Roman', 14),textvariable= d_text2 ,justify='left', width= 9,  bg='white', fg='black')
donor_text2.place(relx=.73, rely=0.70, anchor='w')

#Aceeptors
acceptor_label = Label(filter_frame, text="Acceptors", font=('Times New Roman', 17,'bold'), fg='black', bg=bkg1)
acceptor_label.place(relx=.51, rely=0.75, anchor='e')

a_text1,a_text2 = StringVar(),StringVar()
acceptor_text1 = Entry(filter_frame, font=('Times New Roman', 14),textvariable= a_text1 ,justify='left', width = 9, bg='white', fg='black')
acceptor_text1.place(relx=.52, rely=0.75, anchor='w')

acceptor_text2 = Entry(filter_frame, font=('Times New Roman', 14),textvariable= a_text2 ,justify='left', width= 9,  bg='white', fg='black')
acceptor_text2.place(relx=.73, rely=0.75, anchor='w')

#Rotatable bonds
rbond_label = Label(filter_frame, text="Rotatable Bond Count", font=('Times New Roman', 17,'bold'), fg='black', bg=bkg1)
rbond_label.place(relx=.51, rely=0.80, anchor='e')

rb_text1 , rb_text2 = StringVar(),StringVar()
rbond_text1 = Entry(filter_frame, font=('Times New Roman', 14), textvariable= rb_text1 ,justify='left', width = 9, bg='white', fg='black')
rbond_text1.place(relx=.52, rely=0.80, anchor='w')

rbond_text2 = Entry(filter_frame, font=('Times New Roman', 14),textvariable= rb_text2 ,justify='left', width= 9,  bg='white', fg='black')
rbond_text2.place(relx=.73, rely=0.80, anchor='w')

#Polar_Surface_Area
psa_label = Label(filter_frame, text="Polar Surface Area", font=('Times New Roman', 17,'bold'), fg='black', bg=bkg1)
psa_label.place(relx=.51, rely=0.85, anchor='e')

tpsa_text1, tpsa_text2 = StringVar(),StringVar()
psa_text1 = Entry(filter_frame, font=('Times New Roman', 14),textvariable= tpsa_text1 ,justify='left', width = 9, bg='white', fg='black')
psa_text1.place(relx=.52, rely=0.85, anchor='w')

psa_text2 = Entry(filter_frame, font=('Times New Roman', 14), textvariable= tpsa_text2 ,justify='left', width= 9,  bg='white', fg='black')
psa_text2.place(relx=.73, rely=0.85, anchor='w')

#Fused Aromatic Ring
farc_label = Label(filter_frame, text="Fused Aromatic Rings", font=('Times New Roman', 17,'bold'), fg='black', bg=bkg1)
farc_label.place(relx=.51, rely=0.90, anchor='e')

far_text1,far_text2 = StringVar(),StringVar()
farc_text1 = Entry(filter_frame, font=('Times New Roman', 14),textvariable= far_text1 ,justify='left', width = 9, bg='white', fg='black')
farc_text1.place(relx=.52, rely=0.90, anchor='w')

farc_text2 = Entry(filter_frame, font=('Times New Roman', 14),textvariable= far_text2 ,justify='left', width= 9,  bg='white', fg='black')
farc_text2.place(relx=.73, rely=0.90, anchor='w')


#Add Apply filter Button
apply_filt_button = Button(filter_frame, text = "Customize", font =('Times New Roman', 16,'bold'),bg= bkg , width = 10, command=apply_parameters)
apply_filt_button.place(relx=.58, rely=0.97, anchor='se')

#Add Clear Filter Button
clear_filt_button = Button(filter_frame, text = "Clear", font =('Times New Roman', 16,'bold'),bg= bkg , width = 10, command=resetData)
clear_filt_button.place(relx=.60, rely=0.97, anchor='sw')

# infinite loop 
ws.mainloop()
logger.info(f'GUI Loop Closed at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')

#Close cursor and connection
cur.close()
conn.close()
logger.info(f'Database Connection Closed at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n')

#===============================================================#
