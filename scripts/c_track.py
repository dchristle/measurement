# Now just keep tracking until the process is quit by the user.
track_on = True
fbl_t = qt.instruments['fbl']
track_iter = 0
while track_on == True:
    track_iter = track_iter + 1
    print 'Tracking for %d iteration.' % track_iter
    fbl_t.optimize()
    time.sleep(1.0)
    if msvcrt.kbhit() or track_on == False:
                kb_char=msvcrt.getch()
                if kb_char == "q" or track_on == False: break