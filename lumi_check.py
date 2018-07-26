import fileinput
myfile = open("2011lumibyls.csv","r")
#splitme = [line.split()[0].split(",")[1].split(":") for line in myfile if not line.startswith("#")]
count = 0
for line in myfile:
    if not line.startswith("#"):
        if(line.split(",")[1].split(":")) > 2:
            print "something went wrong", line.split(",")[1]
    count += 1
    if count == 15:
        break


"""
for lsrange in splitme:
    print lsrange, len(lsrange)
    if len(lsrange)>2:
        print "Something went wrong..."
        
"""
