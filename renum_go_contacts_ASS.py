import numpy as np

def _format_line(chain1, atomid1, chain2, atomid2, sigma, epsilon):
    opts = [int(chain1), int(atomid1), int(chain2), int(atomid2), "{:07.5f}".format(sigma), "{:07.5f}".format(epsilon)]
    l = "{0[0]:>10}{0[1]:>10}{0[2]:>10}{0[3]:>10}{0[4]:>15}{0[5]:>15}\n".format(opts)
    return l

arr = np.loadtxt('system.go.parameters_TRI',skiprows=1)
print(arr.shape)
print(len(arr[:,0]))
step = int(22 + 69)
#this is 3384 lines long, you take the length of the first dimension in order to iterate over all of them. 
lines = []
for i in range(len(arr[:,0])):
    chain1   = arr[i,0]
    if chain1 == 1:
        atomid1 = arr[i,1]
    elif chain1 == 2:
        atomid1 = arr[i,1]+step
    elif chain1 == 3:
        atomid1 = arr[i,1]+2*step
    else: 
        print('error')
    chain1 = 1
    chain2   = arr[i,2]
    if chain2 == 1:
        chain2 = 1
        atomid2 = arr[i,3]
    elif chain2 == 2:
        chain2 = 1
        atomid2 = arr[i,3] + step 
    elif chain2 == 3:
        chain2 = 1
        atomid2 = arr[i,3] + 2*step 
    else:
        print('error2')
    chain2 = 1
    sigma    = arr[i,4]
    epsilon  = arr[i,5]
    lines.append(_format_line(chain1,atomid1,chain2,atomid2,sigma,epsilon))

with open('system.go.parameters_TRI_converted', 'w+') as outfile:
    outfile.write("12-10 potentials\n")
    for line in lines:
        outfile.write(line)
