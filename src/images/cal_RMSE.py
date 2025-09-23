input_num,mutation_num = map(int,input().split())

true_j = []
est_j = []
for i in range(input_num):
    true_j.append(float(input()))
input()
for i in range(input_num):
    est_j.append(float(input()))

R = 0
cnt = 0
for i in range(input_num):
    #if true_j[i]>0.95:
        R += (true_j[i] - est_j[i])**2
        cnt += 1
    #print(abs(true_j[i] - est_j[i]))
    #if abs(true_j[i] - est_j[i])>0.04:
        #print(true_j[i], est_j[i])

R /= cnt
print(R**0.5)