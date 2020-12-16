def files(data_folder):

    import os 

    os.chdir(data_folder)
    folder_list=[]
    file_list=[]
    for i in os.listdir():
        if i.startswith("._")==False:
            folder_list.append(i)

    for j in folder_list:
        os.chdir(data_folder + "/" + j)
        for k in os.listdir():
            if k.startswith("._")==False:
                file_list.append(j+ "/" +k)

                # with open(data_folder + "/" + "DataFiles.txt","a", encoding="utf-8") as file:
                #     file.write(k+"\n")

    return file_list

# files("/Volumes/Elements/Dersler/Master Courses/EQE 520/Assignment 1/Python/data")

    


    # for j in i:
    #     with open(j, "r", encoding="utf-8") as file:
    #         test_value="STATION CODE: "

    #         file.seek(file.read().find(test_value)+len(test_value))
    #         interval=file.read(4)
    #         interval=[''.join(interval)]
    #         dt=interval[0]
    #         print(dt)

