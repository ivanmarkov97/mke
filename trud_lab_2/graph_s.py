import matplotlib.pyplot as plt

data_file = open("data_s.txt")
y = []
x = []
graph_type = data_file.readline()
str_num = data_file.readline()
num_el = int(str_num)
step = 45.0 / (num_el)

for i in range(num_el + 1):
	y.append(step*i)
for line in data_file:
	line = line.replace('\n', '')
	line = line.replace(',', '')
	x.append(float(line))

print(y)
print(x)

plt.plot(y,x)
plt.xlabel("x")
plt.ylabel("function(x)")
plt.text(2, 3, "diff == " + str( (x[-1] - x[-2]) / (y[-1] - y[-2])))
plt.text(2, 2.8, "elem size == " + str(step))
plt.show()

