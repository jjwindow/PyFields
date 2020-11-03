from all_funcs import field_trace, dipole
import matplotlib.pyplot as plt

pstart_1 = [1., 0.1, 0.1]
pstart_2 = [1., 0.15, 0.1]
pstart_3 = [1., 0.2, 0.1]
pstart_4 = [1., 0.25, 0.1]
# p_dip_1, B_dip_1 = field_trace(pstart_1, dipole, 0.01, 100000)
# p_dip_2, B_dip_2 = field_trace(pstart_2, dipole, 0.01, 100000)

# plt.polar(p_dip_1[:, 1], p_dip_1[:, 0])
# plt.polar(p_dip_2[:, 1], p_dip_2[:, 0])

x_1, y_1 = field_trace(pstart_1, dipole, 0.01, 100000)
# x_2, y_2 = field_trace(pstart_2, dipole, 0.01, 100000)
# x_3, y_3 = field_trace(pstart_3, dipole, 0.01, 100000)
# x_4, y_4 = field_trace(pstart_4, dipole, 0.01, 100000)

plt.xlabel("x (arbitrary units)")
plt.ylabel("y (arbitrary units")
plt.plot(x_1, y_1, color = 'b')
# plt.plot(x_2, y_2, color = 'b')
# plt.plot(x_3, y_3, color = 'b')
# plt.plot(x_4, y_4, color = 'b')
plt.title("Test Dipole")
plt.show()
#print(f"R:\n{r[:20]}")
#print(f"THETA:\n{theta[:20]}")