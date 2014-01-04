import datetime


def get_timestamp():
    return datetime.datetime.now().strftime("%Y%m%d-%H%M%S")


def get_all_module_names(dirname):
    dirname = os.path.abspath(dirname)
    names = []
    for root, _, files in os.walk(dirname):
        for file in files:
            if not file.endswith('.py') or file == '__init__.py' or file == 'base.py':
                continue
            name = os.path.join(root, file)[len(dirname)+1:].replace('/','.')
            names.append(name[:-3])
    return sorted(names)
