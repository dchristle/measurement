import time
import progressbar

# with progressbar.ProgressBar(max_value=10) as bar:
#     for i in range(10):
#         time.sleep(0.1)
#         bar.update(i)

bar = progressbar.ProgressBar(redirect_stdout=False)
for i in range(100):
    print 'Some text', i
    time.sleep(0.1)
    bar.update(i)