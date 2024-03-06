import sys
import time

def _print_transcript_real_time(filename, frequency=0.1):
    """Stream new lines from a file in real time."""
    with open(filename, 'r') as file:
        # Move the pointer to the end of the file
        file.seek(0, 2)
        while True:
            line = file.readline()
            if line:
                print('\t', line, end='')
            else:
                time.sleep(frequency)

if __name__ == "__main__":

    # Get the transcript file as command-line argument
    if len(sys.argv) != 2:
        print("Usage: python stream_output.py <file_path>")
        sys.exit(1)

    # Display transcript as the file gets updated  
    _print_transcript_real_time(filename=sys.argv[1])
