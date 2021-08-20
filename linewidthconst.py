

class MetaChoise(type):
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(MetaChoise, cls).__call__(*args, **kwargs)
        return cls._instances[cls]

class LineWidthConst(metaclass = MetaChoise):
    
    pass

if __name__ == '__main__':
    s = LineWidthConst()
    print(id(s))
    m = LineWidthConst()
    print(id(m))
    print(id(s)==id(m))

