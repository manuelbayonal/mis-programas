class triangulo:
    def __init__ (self,lado,altura):
        self.lado=lado
        self.altura=altura
    def area(self):
        return (self.lado*self.altura)/2 
if __name__=='__main__':
    triangulo=triangulo(lado=2,altura=4)
    print("el área del triangulo es: "+str(triangulo.area()))