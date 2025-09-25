import sys
import os
import pypulseq as pp

# Predefined acoustic resonance frequencies (example values - adjust as needed)
DEFAULT_ACOUSTIC_RESONANCES = [ # For our Siemens 3T Prisma
    {"frequency": 590, "bandwidth": 100},
    {"frequency": 1140, "bandwidth": 220}
]

def check_PNS(seq, gradFile):
    """
    Plots the PNS level of the sequence using the gradient information in gradFile and throws an error if PNS level is exceeded
    """
    import matplotlib.pyplot as plt

    # Check if gradient file exists
    gradient_file = os.path.join(os.path.dirname(__file__), gradFile)
    if os.path.exists(gradient_file):
        print("Gradient file found, using it for PNS check.")
    else:
        print("Gradient file not found, skipping PNS check.")
        return

    pns_ok, pns_n, pns_c, tpns = seq.calculate_pns(gradient_file, do_plots=True)
    if pns_ok:
        print("PNS check passed successfully: maximum PNS level =", 100*pns_n.max(), "%")
        plt.show()
    else:
        print("PNS check failed: maximum PNS level =", 100*pns_n.max(), "%")
        plt.show()
        raise RuntimeError("PNS check failed!")

def check_acoustic_resonances(seq, gradFile):
    """
    Plots the acoustic resonances of the sequence against the forbidden frequencies specified in gradFile
    """
    import matplotlib.pyplot as plt
    gradient_file = os.path.join(os.path.dirname(__file__), gradFile) # Path to the Prisma gradient file
    
    if os.path.exists(gradient_file):
        # Use gradient file if available
        from pypulseq.utils.siemens.readasc import readasc
        asc, extra = readasc(gradient_file) # Parse the ASC file

        # Extract acoustic resonances
        frequencies = asc["asGPAParameters"][0]["sGCParameters"]["aflAcousticResonanceFrequency"]
        bandwidths  = asc["asGPAParameters"][0]["sGCParameters"]["aflAcousticResonanceBandwidth"]
        #print("Frequencies:", frequencies)
        #print("Bandwidths:", bandwidths)

        # Combine frequencies and bandwidths into a list of dictionaries
        acoustic_resonances = [
            {"frequency": freq, "bandwidth": bw}
            for freq, bw in zip(frequencies.values(), bandwidths.values()) # Use .values() to get the dictionary values
            if freq != 0  # Exclude entries with zero frequency
        ]
        print(f"Using acoustic resonances from gradient file: {gradFile}")
    else:
        # Use default acoustic resonances if gradient file is not available
        acoustic_resonances = DEFAULT_ACOUSTIC_RESONANCES
        print(f"Warning: Gradient file '{gradFile}' not found. Using default acoustic resonances.")

    seq.calculate_gradient_spectrum(acoustic_resonances=acoustic_resonances, use_derivative=True, frequency_oversampling=10, window_width=0.5)
    plt.title("Acoustic Resonances")
    plt.show()

def safety_check(seq_file):
    """
    Test a pulseq sequence for PNS and Acoustic Resonances
    """
    gradFile = "MP_GPA_K2309_2250V_951A_AS82.asc" # Siemens 3T Prisma gradient file for PNS and acoustic resonance check (place in the same folder as this file)

    # Create a Sequence object
    seq = pp.Sequence()

    # Load the sequence from a .seq file
    seq.read(seq_file)

    # Check the PNS level
    check_PNS(seq, gradFile)

    # Check the acoustic resonances
    check_acoustic_resonances(seq, gradFile)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python safety_check.py <path_to_seq_file>")
        sys.exit(1)

    seq_file = sys.argv[1]
    safety_check(seq_file)