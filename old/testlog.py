import logging
import logging.handlers
import multiprocessing
from time import sleep
from random import random, randint


# Almost the same as the demo code, but added `console_handler` to directly
# read logging info from the console
def listener_configurer():
    root = logging.getLogger()
    file_handler = logging.handlers.RotatingFileHandler('mptest.log', 'a', 300, 10)
    console_handler = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s %(processName)-10s %(name)s %(levelname)-8s %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    root.addHandler(file_handler)
    root.addHandler(console_handler)
    root.setLevel(logging.DEBUG)


# Almost the same as the demo code, but made it into a forever loop. This is
# more likely to happen in an app that does not have a clear end point, e.g.
# a deployed IoT sensor. Another change is to show that if configurer is not
# passed but directly visible by the process function, calling it directly has
# the same effect.
def listener_process(queue):
    listener_configurer()
    while True:
        while not queue.empty():
            record = queue.get()
            logger = logging.getLogger(record.name)
            logger.handle(record)  # No level or filter logic applied - just do it!
        sleep(1)


# Same as demo code
def worker_configurer(queue):
    h = logging.handlers.QueueHandler(queue)  # Just the one handler needed
    root = logging.getLogger()
    root.addHandler(h)
    # send all messages, for demo; no other level or filter logic applied.
    root.setLevel(logging.DEBUG)


# Almost the same as demo code, except the logging is simplified, and configurer
# is no longer passed as argument.
def worker_process(queue):
    worker_configurer(queue)
    for i in range(3):
        sleep(random())
        innerlogger = logging.getLogger('worker')
        innerlogger.info(f'Logging a random number {randint(0, 10)}')


def main():
    queue = multiprocessing.Queue(-1)
    listener = multiprocessing.Process(
        target=listener_process, args=(queue,))
    listener.start()
    workers = []
    for i in range(3):
        worker = multiprocessing.Process(target=worker_process, args=(queue,))
        workers.append(worker)
        worker.start()
    for w in workers:
        w.join()


if __name__ == '__main__':
    main()
