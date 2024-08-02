band_freq_range_dict = {
"UWB1": [500, 650, 800, 950],
"UWB2": [950, 1150, 1350, 1550, 1750],
"UWB3": [1750, 1950, 2150, 2350, 2550],
"UWB4": [2550, 2700, 2850, 3000, 3150, 3300, 3450]
}

def get_band_range(rest_freq: float):
    """
    get the corresponding band and frequency range according the given rest frequency.
    
    Parameter
    ==========
    rest_freq: float

    Return
    ======
    return: (str, [int, int])
    """
    for band_name, freq_intervals in band_freq_range_dict.items():
        for i in range(len(freq_intervals)-1):
            start_freq, end_freq = freq_intervals[i], freq_intervals[i+1]
            if start_freq <= rest_freq <= end_freq:
                return (band_name, [start_freq, end_freq])
    return ("", [0, 0])

if __name__ == "__main__":
    band_name, freq_range = get_band_range(1665)
    print(band_name, freq_range)

