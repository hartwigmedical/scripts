def perform_prod_test(profile):
    if profile == 'prod':
        prod_warn = input("Warning: you are running in prod. Type 'y' to continue.\n")
        if prod_warn.lower() != 'y':
            print('Program aborted')
            exit(1)