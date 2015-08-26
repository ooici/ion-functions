def version(version_number):
    """Wrapper to add version to ion function"""
    def decorate(f):
        f.version = version_number
        return f
    return decorate