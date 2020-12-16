def nonVs30(path):
    import datatolist
    import os

    # path="/Volumes/Elements/Dersler/Master Courses/EQE 520/Assignment 1/Python/data"
    files=datatolist.files(path)
    for i in range(0, len(files)):
        with open(path + "/" + files[i], "r", encoding="utf-8") as file:
            test_value="VS30_M/S: "

            file.seek(file.read().find(test_value)+len(test_value))
            Vs30_value=file.read(4)
            Vs30_value=[''.join(Vs30_value)]
            Vs30_value=Vs30_value[0]
            if Vs30_value=="None":
                os.remove(path + "/" + files[i])

def streamCode(path):
    import datatolist
    import os
    
    with open(path, "r", encoding="utf-8") as file:
    
        test_value2="STREAM: "

        file.seek(file.read().find(test_value2)+len(test_value2))
        stream_code=file.read(3)
        stream_code=[''.join(stream_code)]
        stream_code=stream_code[0]
    
    return stream_code

