p = 998244353
g = 3
g_inv = pow(g, p-2, p)
max_n = 2048

# File path in notebook environment
file_path = './roots.txt'
output_file = './inv_roots.txt'

with open(file_path, 'w') as f:
    # Line 0: dummy for n=0 (unused)
    f.write('00000000\n')
    # Lines 1..1000: principal n-th root of unity or 0 if not supported
    for n in range(1, max_n + 1):
        if (p - 1) % n == 0:
            omega = pow(g, (p - 1) // n, p)
        else:
            omega = 0
        # Write as 8-digit hex (zero-padded)
        f.write(f'{omega:08x}\n')

print(f"Wrote principal roots into: {file_path}")

with open(output_file, 'w') as f:
    # Line 0: placeholder for n=0 (unused)
    f.write('00000000\n')
    for n in range(1, max_n + 1):
        if (p - 1) % n == 0:
            # compute inverse root: (g_inv)^{(p-1)/n} mod p
            inv_omega = pow(g_inv, (p - 1) // n, p)
        else:
            inv_omega = 0
        # write as 8-digit zero-padded hex
        f.write(f'{inv_omega:08x}\n')

print(f"Wrote inverse roots up to n={max_n} into '{output_file}'")

output_file = 'inv_n.txt'

with open(output_file, 'w') as f:
    # Line 0: placeholder for n=0 (unused)
    f.write('00000000\n')
    for n in range(1, max_n + 1):
        inv_n = pow(n, p-2, p)
        f.write(f'{inv_n:08x}\n')

print(f"Wrote modular inverses 1..{max_n} to '{output_file}'")
