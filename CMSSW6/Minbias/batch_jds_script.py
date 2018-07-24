def makeBatchConfigjdsFile(job_dir):

    config  = 'executable = run.sh \n'
    config += 'universe   = vanilla \n'
    config += 'log        = job_$(Process).log \n'
    config += 'output     = job_$(Process).out \n'
    config += 'error      = job_$(Process).err \n'
    config += 'transfer_input_files = x_test.C \n'
    config += 'should_transfer_files = YES  \n'
    config += 'when_to_transfer_output = ON_EXIT  \n'
    config += 'queue'

    return config











