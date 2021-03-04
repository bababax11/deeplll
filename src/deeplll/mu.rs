use rug::Rational;
#[derive(PartialEq, Eq, Clone, Debug)]
pub struct Mu {
    arr: Vec<Rational>,
}

impl std::ops::Index<(usize, usize)> for Mu {
    type Output = Rational;
    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let (a, b) = index;
        &self.arr[a * (a - 1) / 2 + b]
    }
}
impl std::ops::IndexMut<(usize, usize)> for Mu {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let (a, b) = index;
        &mut self.arr[a * (a - 1) / 2 + b]
    }
}

impl Mu {
    pub fn new(n: usize) -> Self {
        Mu {
            arr: vec![Rational::from(0); n * (n - 1) / 2],
        }
    }
}
