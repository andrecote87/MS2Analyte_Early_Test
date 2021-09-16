import fileinput
import sys
import os









# string1 = 'intensity_cutoff = '
  
# # opening a text file
# file1 = open("config.py", "r")
  
# # setting flag and index to 0
# flag = 0
# index = 0
  
# # Loop through the file line by line
# for line in file1:  
#     index += 1 
#     print(line)
#     # checking string is present in line or not
#     if string1 in line:
      
#       flag = 1
#       break 

# # checking condition for string found or not
# if flag == 0: 
#    print('String', string1 , 'Not Found') 
# else: 
#    print('String', string1, 'Found In Line', index)
#    print(line)

  
# # closing text file    
# file1.close() 




# # opening the file in read mode
# file = open("config_test.py", "r")
# all_lines = file.readlines()
# print(all_lines[78])
# replacement = ""
# # using the for loop
# for line in file:
#     line = line.strip()
#     changes = line.replace(all_lines[78], "intensity_cutoff = 3000")
#     replacement = replacement + changes + "\n"

# file.close()
# # opening the file in write mode
# fout = open("config_test.py", "w")
# fout.write(replacement)
# fout.close()


absolute_path = os.path.abspath(__file__)
print("Directory Path: ", os.path.dirname(absolute_path))
file2 = open(os.path.join(os.path.dirname(absolute_path),'config_test.py'), "r")

file = open('config_test.py')
all_lines = file2.readlines()
print(all_lines[78])


all_lines[78] = "intensity_cutoff = 8000\n"





with open("config_test.py", 'w') as output:
    for row in all_lines:
        output.write(str(row))







# absolute_path = os.path.abspath(__file__)
# print("Full path: ", absolute_path)
# print("Directory Path: ", os.path.dirname(absolute_path))
# print(self.ThresholdLineEdit.text())





# # file = open(os.path.dirname(absolute_path)+"\\config.py")

# # print(file)
# # # print('test')
# # all_lines = file.readlines()
# # print(all_lines[78])




# # for line in fileinput.input(file, inplace=1):
# #     print(line)
# #     line = line.replace(all_lines[78], "intensity_cutoff = " + self.ThresholdLineEdit.text() +'\n')
# #     sys.stdout.write(line)
# # # var1 = all_lines[78]
# # # var2 = "intensity_cutoff = " + self.ThresholdLineEdit.text() +'\n'
# # # file = "config.py"
# # # replacement(file, var1, var2)




# # opening the file in read mode
# file1 = open(os.path.join(os.path.dirname(absolute_path),'config.py'), "r")
# print(file1)
# file2 = open('config.py','r')
# print(file2)

# # print('test')
# all_lines = file2.readlines()
# print('test')
# replacement = ""
# print('test')
# # using the for loop
# for line in file2:
#     print('test')
#     line = line.strip()
#     print('test')
#     changes = line.replace(all_lines[78], "intensity_cutoff = 3000")
#     print('test')
#     replacement = replacement + changes + "\n"

# file2.close()
# # opening the file in write mode
# fout = open("config_test.py", "w")
# fout.write(replacement)
# fout.close()