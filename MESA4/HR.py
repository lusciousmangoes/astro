

#change this for each folder
n=66 #maximum profile number in folder

data = open('HR.dat','a')

for i in range(1,n+1): 
    f = open("profile"+str(i)+".data","r")
    lines = f.readlines()

	# If we need to read line 33, and assign it to some variable
    x = lines[6]
    h=x.split()
    nums = [h[2],h[5]]
    for n in range(0,len(nums)):
        nums[n] = float(nums[n])
        if n%2==0:
            print >> data, nums[n], str(10.**float(nums[n-1]))
            # luminosity, temperature

data.close()
