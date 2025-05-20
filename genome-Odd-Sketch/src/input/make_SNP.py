import random

num_input,num_mutation = map(int,input().split())

print(64,num_input)

ATGC = ["A","T","G","C"]
inputs = [[] for _ in range(num_input)]

for i in range(num_input):
    input1 = []
    input2 = []
    len_genome = 5*(10**5) + random.randint(-10**5,10**5)
    
    for j in range(len_genome):
        input1.append(ATGC[random.randint(0,3)])
        input2.append(input1[-1])
    
    mut = random.randint(10,num_mutation)
    p = random.sample(range(len_genome-1), k=mut)
    
    for j in p:
        input2[j] = ATGC[(random.randint(1,3)+ATGC.index(input2[j]))%4] 
        
    inputs[i].append("".join(input1))
    inputs[i].append("".join(input2))
    print(len_genome)

for i in range(num_input):
    print(inputs[i][0])
    print(inputs[i][1])