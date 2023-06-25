void test(){

struct compare
{
    int key;
    compare(int const &i): key(i) {}
 
    bool operator()(int const &i) {
        return (i == key);
    }
};

 std::vector<int> v = { 7, 3, 6, 2, 6 };

 int run = 6; // user input run
 int run_idx;
 
 // search for run number
 auto itr = std::find_if(v.cbegin(), v.cend(), compare(run));

 // if found run, get the index
 if (itr != v.cend()) {
   run_idx = std::distance(v.cbegin(), itr);
   cout << "run index:" << run_idx << endl;
   //std::cout << "Element present at index " << std::distance(v.cbegin(), itr) << endl;
 }
 else {
   std::cout << "run number not found" << endl;
   
 }
 
}
