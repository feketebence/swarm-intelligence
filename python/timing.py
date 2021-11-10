import time
import datetime

def echo_running_time():
    current_timestamp = time.time()
    date_time = datetime.datetime.fromtimestamp(current_timestamp)

    print(f"Run finished on: ", date_time)

def main():
    echo_running_time()

if __name__ == "__main__":
    main()
