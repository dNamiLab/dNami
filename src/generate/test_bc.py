from genRhs import incPATH
bc_template = open(incPATH[:-13]+'template_genbc.for','r')

input  = open(incPATH[:-13]+'rhs.for','r')
output = open(incPATH[:-13]+'rhs.for','w')

lines  = input.readlines()


# bctype = 10

# edit = [[],[]]
# edit[0].append('!PHYBCLOCVAR');edit[1].append('\# include include_PhyBClocVar_'+bcname+'_'+dirBC+'_'+bctype+'.f90')

# for bctype in range(10,13):
# 	for l in lines:
# 		outpout.write(l.replace('BCNAME','bc_type'+str(bctype)))

# 	outpout.write("\n\n")


for l in lines:
	output.write(l.replace('cmpstored','testreplace'))
