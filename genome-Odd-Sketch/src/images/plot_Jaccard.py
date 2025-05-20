import matplotlib.pyplot as plt
import numpy as np
input_num,mutation = map(int,input().split())
true_J = []
cal_J = []
for i in range(input_num):
    true_J.append(float(input()))
input()
for i in range(input_num):
    cal_J.append(float(input()))
    
x = np.linspace(0, 1, 100)
y = x

plt.xlim(0.60,1.00)
plt.ylim(0.60,1.00)
plt.xlabel('Jaccard_true')
plt.ylabel('Jaccard_estimate')
c1, c2 = 'blue', 'red'

plt.scatter(true_J,cal_J, color=c1, zorder=2)
plt.plot(x,y, color=c2, zorder=1)


plt.show()
    