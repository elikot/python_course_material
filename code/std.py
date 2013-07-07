def std(foo, verbose=False): # foo is a list
    total=0
    for n in foo:
        total+=n
    average=float(total)/len(foo) # nb! float, et poleks int vastuseks
    A=0
    for item in foo:
        A+=(item-average)**2
        #print(A)
        stdev=(A/len(foo))**0.5
    if verbose:
        print('Total: {0}.\n'.format(total))
        print('Average: {0}.\n'.format(average))
        print('Std: {0}.\n'.format(stdev))
    return stdev
