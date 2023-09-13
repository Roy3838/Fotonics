# Given wavelengths in nm
l1 = 871
l2 = 980
l3 = 1064

# Calculate possible generated wavelengths
wavelengths = {
    "2l1": 2 * l1,
    "2l2": 2 * l2,
    "2l3": 2 * l3,
    "l1 + l2": l1 + l2,
    "l1 - l2": abs(l1 - l2),
    "l1 + l3": l1 + l3,
    "l1 - l3": abs(l1 - l3),
    "l2 + l3": l2 + l3,
    "l2 - l3": abs(l2 - l3)
}

print(wavelengths)