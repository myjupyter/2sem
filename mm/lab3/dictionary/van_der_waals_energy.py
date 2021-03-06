from dictionary.process import process_van_der_waals_energy

__van_der_waals_energy = {
    "data": {
        "H ": {"R*j": 0.6000, "e_k": 0.0157},
        "HO": {"R*j": 0.0000, "e_k": 0.0000},
        "HS": {"R*j": 0.6000, "e_k": 0.0157},
        "HC": {"R*j": 1.4870, "e_k": 0.0157},
        "H1": {"R*j": 1.3870, "e_k": 0.0157},
        "H2": {"R*j": 1.2870, "e_k": 0.0157},
        "H3": {"R*j": 1.1870, "e_k": 0.0157},
        "HP": {"R*j": 1.1000, "e_k": 0.0157},
        "HA": {"R*j": 1.4590, "e_k": 0.0150},
        "H4": {"R*j": 1.4090, "e_k": 0.0150},
        "H5": {"R*j": 1.3590, "e_k": 0.0150},
        "HW": {"R*j": 0.0000, "e_k": 0.0000},
        "O ": {"R*j": 1.6612, "e_k": 0.2100},
        "O2": {"R*j": 1.6612, "e_k": 0.2100},
        "OW": {"R*j": 1.7682, "e_k": 0.1521},
        "OH": {"R*j": 1.7210, "e_k": 0.2104},
        "OS": {"R*j": 1.6837, "e_k": 0.1700},
        "CT": {"R*j": 1.9080, "e_k": 0.1094},
        "C ": {"R*j": 1.9080, "e_k": 0.0860},
        "N ": {"R*j": 1.8240, "e_k": 0.1700},
        "N3": {"R*j": 1.8240, "e_k": 0.1700},
        "S ": {"R*j": 2.0000, "e_k": 0.2500},
        "SH": {"R*j": 2.0000, "e_k": 0.2500},
        "P ": {"R*j": 2.1000, "e_k": 0.2000},
        "IM": {"R*j": 2.47, "e_k": 0.1},
        "Li": {"R*j": 1.1370, "e_k": 0.0183},
        "IP": {"R*j": 1.8680, "e_k": 0.00277},
        "K ": {"R*j": 2.6580, "e_k": 0.000328},
        "Rb": {"R*j": 2.9560, "e_k": 0.00017},
        "Cs": {"R*j": 3.3950, "e_k": 0.0000806},
        "I ": {"R*j": 2.35, "e_k": 0.40},
        "F ": {"R*j": 1.75, "e_k": 0.061},
        "IB": {"R*j": 5.0, "e_k": 0.1},
        "EP": {"R*j": 0.0, "e_k": 0.0},
    }
}
van_der_waals_energy = process_van_der_waals_energy(__van_der_waals_energy)
