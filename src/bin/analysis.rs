use deeplll::parse::cat;
use once_cell::sync::Lazy;
use regex::Regex;

fn parse_sec_len_cnt(s: &str) -> (f32, u32, u32) {
    static RE_SEC: Lazy<Regex> = Lazy::new(|| Regex::new(r"(.[0-9])+ sec").unwrap());
    static RE_LEN_CNT: Lazy<Regex> =
        Lazy::new(|| Regex::new(r"\(hist.len, cnt\): \([0-9]+, [0-9]+\)").unwrap());
    let sec = &RE_SEC.captures(s).unwrap()[0];
    dbg!(sec);
    let mut sec = sec.split(' ');
    let _ = sec.next();
    let sec = sec.next().unwrap().parse::<f32>().unwrap();

    let len_cnt = &RE_LEN_CNT.captures(s).unwrap()[0];

    static RE_TUP: Lazy<Regex> = Lazy::new(|| Regex::new(r"[0-9]+").unwrap());
    let len_cnt = &mut RE_TUP.captures_iter(len_cnt);
    let (len, cnt) = (&len_cnt.next().unwrap()[0], &len_cnt.next().unwrap()[0]);
    (
        sec,
        len.parse::<u32>().unwrap(),
        cnt.parse::<u32>().unwrap(),
    )
}

fn main() {}

#[cfg(test)]
mod tests {
    use super::*;
    use average::assert_almost_eq;
    #[test]
    fn parse_sec_len_cnt_test() {
        let s = cat("results/deeplll/ndim10seed0cnt1000dim10delta0.98999995.txt".as_ref()).unwrap();
        let (sec, len, cnt) = parse_sec_len_cnt(&s);

        assert_almost_eq!(sec, 3.951, 1e-5);
        assert_eq!(len, 486);
        assert_eq!(cnt, 6316);
    }
}