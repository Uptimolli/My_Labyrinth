#include <iostream>
#include <stack>
#include <vector>
#include <string>
#include <queue>

using namespace std;


class Bool_Matrix {
private:
    int n; int m;

    vector<vector<bool>> my_matrix;
public:
    
    Bool_Matrix (int g_n = 3, int g_m = 3) {
        vector<vector<bool>>  new_matrix(g_n, vector <bool>(g_m));
        for (int i = 0; i < g_n; i++)
            for (int j = 0; j < g_m; j++)
                new_matrix[i][j] = 0;
        my_matrix = new_matrix;
        n = g_n;
        m = g_m;
        
    }
    int get_n() {
        return n;
    }
    

    int get_m() {
        return m;
    }

    bool get_kl(int k, int l) {
            return my_matrix[k][l];
    }

    void set_new_value(int index_x, int index_y, bool new_value) {
        my_matrix[index_x][index_y] = new_value;
    }


    void output_k_l_element(int k, int l) {
        cout<< my_matrix[k][l]<< ' ';
    }
    void output_k_stroka(int k) {
        for (int i = 0; i < m; i++)
            cout << my_matrix[k][i]<<' ';
        cout << endl;
    }
    void output_matrix() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++)
                cout << my_matrix[i][j] << ' ';
            cout << endl;
        }
    }
};


class Labyrinth {
private:
    vector<string> potw{ "UP", "DOWN", "LEFT", "RIGHT" };
    vector<int> start;
    vector<int> exit;
    Bool_Matrix laby;

    int big_number = 1000000000;

    int do_a_correct_number(int k) {
        return 2 * k + 1;
    }
    int do_a_real_number(int k) {
        return (k - 1) / 2;
    }
public:
    Labyrinth(int n, int m) {
        n = do_a_correct_number(n);
        m = do_a_correct_number(m);

        Bool_Matrix lab(n, m);
        laby = lab;
    }

    void create_2d_lab() {
        srand(time(NULL));
        try {
            int laby_n = laby.get_n();
            int laby_m = laby.get_m();
            // Делаем случайную точку началом для крота и делаем координаты нечетными
            int actualX = rand() % laby_n; actualX = actualX - (actualX + 1) % 2;
            int actualY = rand() % laby_m; actualY = actualY - (actualY + 1) % 2;
            if (actualX == -1) actualX = 1; if (actualY == -1) actualY = 1;
            laby.set_new_value(actualX, actualY, 1);
            

            string where_we_would_go;

            int unvisit = do_a_real_number(laby_n)*do_a_real_number(laby_m) - 1;

            vector<vector<bool>> visited(laby_n, vector<bool> (laby_m));


            visited[actualX][actualY] = true;
            

            while (unvisit) {
                where_we_would_go = potw[rand() % 4];
                if (where_we_would_go == "LEFT" && actualY > 1) {
                    if (!visited[actualX][actualY - 2]) {
                        visited[actualX][actualY - 1] = true;
                        visited[actualX][actualY - 2] = true;
                        unvisit--;
                    }
                    if (visited[actualX][actualY - 1] &&
                        visited[actualX][actualY - 2])
                    {
                        actualY = actualY - 2;
                    }
                }
                if (where_we_would_go == "RIGHT" && actualY < laby_m - 2) {
                    if ( !visited[actualX][actualY + 2]) {
                        visited[actualX][actualY + 1] = true;
                        visited[actualX][actualY + 2] = true;
                        unvisit--;
                    }
                    if (visited[actualX][actualY + 1] &&
                        visited[actualX][actualY + 2])
                    {
                        actualY = actualY + 2;
                    }
                }
                if (where_we_would_go == "UP" && actualX > 1) {
                    if (!visited[actualX - 2][actualY]) {
                        visited[actualX - 1][actualY] = true;
                        visited[actualX - 2][actualY] = true;
                        unvisit--;
                    }
                    if (visited[actualX - 1][actualY] &&
                        visited[actualX - 2][actualY])
                    {
                        actualX = actualX - 2;
                    }

                }
                if (where_we_would_go == "DOWN" && actualX < laby_n - 2) {
                    if ( !visited[actualX + 2][actualY]) {
                        visited[actualX + 1][actualY] = true;
                        visited[actualX + 2][actualY] = true;
                        unvisit--;
                    }
                    if (visited[actualX + 1][actualY] &&
                        visited[actualX + 2][actualY])
                    {
                        actualX = actualX + 2;
                    }
                }

            }
            for (int i = 0; i < laby_n; i++)
                for (int j = 0; j < laby_m; j++)
                    laby.set_new_value(i, j, visited[i][j]);
            create_a_random_wall();
            create_a_random_hole();
            create_a_random_hole();
            create_a_random_hole();
        }
        catch(...){
            cout << "Error at creating 2d labyrinth";
        }
    }
    

    void output_laby() {
        //laby.output_matrix();
        for (int i = 0; i < laby.get_n() ; i ++) {
            for (int j = 0; j < laby.get_m(); j ++) {
                laby.output_k_l_element(i, j);
            }
            cout << endl;
        }
    }
    
    
    void create_entry_and_exit() {

        bool entry_flag = true; bool exit_flag = true;
        int exact_x_entry_coordinate, exact_y_entry_coordinate, exact_x_exit_coordinate, exact_y_exit_coordinate;

        while (entry_flag || exit_flag) {
            string correct_place = potw[rand() % 4];
            int exact_x_coordinate, exact_y_coordinate;
            if (correct_place == "UP") {
                exact_x_coordinate = 1 + (rand() % (laby.get_n() - 2));
                exact_x_coordinate = exact_x_coordinate - (exact_x_coordinate + 1) % 2;
                exact_y_coordinate = 0;
            }
            if (correct_place == "DOWN") {
                exact_x_coordinate = 1 + (rand() % (laby.get_n() - 2));
                exact_x_coordinate = exact_x_coordinate - (exact_x_coordinate + 1) % 2;
                exact_y_coordinate = laby.get_m() - 1;
            }

            if (correct_place == "LEFT") {
                exact_x_coordinate = 0;
                exact_y_coordinate = 1 + (rand() % (laby.get_m() - 2));
                exact_y_coordinate = exact_y_coordinate - (exact_y_coordinate + 1) % 2;
            }
            if (correct_place == "RIGHT") {
                exact_x_coordinate = laby.get_n() - 1;
                exact_y_coordinate = 1 + (rand() % (laby.get_m() - 2));
                exact_y_coordinate = exact_y_coordinate - (exact_y_coordinate + 1) % 2;
            }
            if (entry_flag) {
                exact_x_entry_coordinate = exact_x_coordinate;
                exact_y_entry_coordinate = exact_y_coordinate;
                entry_flag = false;
            }
            if (!entry_flag && exact_x_coordinate != exact_x_entry_coordinate && exact_y_entry_coordinate != exact_y_coordinate) {
                exact_x_exit_coordinate = exact_x_coordinate;
                exact_y_exit_coordinate = exact_y_coordinate;
                exit_flag = false;
            }
        }
        start = { exact_x_entry_coordinate, exact_y_entry_coordinate };
        exit = { exact_x_exit_coordinate, exact_y_exit_coordinate };
        laby.set_new_value(exact_x_entry_coordinate, exact_y_entry_coordinate, true);
        laby.set_new_value(exact_x_exit_coordinate, exact_y_exit_coordinate, true);
    }
    
    void insert_something_between_two_points_2d(int x_1, int y_1, int x_2, int y_2, bool something) {
        // 1 is a hole, 0 is a wall
        if (abs(x_1 - x_2) + abs(y_1 - y_2) > 1) {
            cout<<"Ошибка ввода, вы выбрали не рядомстоящие точки"<<endl;
        }
        else {
            int X;
            int Y;
            X = (do_a_correct_number(x_1)+ do_a_correct_number(x_2))/2;
            Y = (do_a_correct_number(y_1)+ do_a_correct_number(y_2))/2;
            laby.set_new_value(X, Y, something);
        }

    }

    void create_a_random_hole() {
        int randomX; int randomY; int secondX; int secondY;
        
        randomX = rand() % do_a_real_number(laby.get_n());
        randomY = rand() % do_a_real_number(laby.get_m());
        
        secondX = randomX;
        secondY = randomY;

        while (!(abs(secondX-randomX+secondY-randomY))){
            string where_we_would_do_a_hole = potw[rand() % 4];
            if (where_we_would_do_a_hole == "UP" && randomX != 0) {
                secondX--;
            }
            if (where_we_would_do_a_hole == "DOWN" && randomX != do_a_real_number(laby.get_n())-1) {
                secondX++;
            }
            if (where_we_would_do_a_hole == "LEFT" && randomY != 0) {
                secondY--;
            }
            if (where_we_would_do_a_hole == "RIGHT" && randomY != do_a_real_number(laby.get_m())-1) {
                secondY++;
            }
        }
        insert_something_between_two_points_2d(randomX, randomY, secondX, secondY, true);

    }
    void create_a_random_wall() {
        int randomX; int randomY; int secondX; int secondY;

        randomX = rand() % do_a_real_number(laby.get_n());
        randomY = rand() % do_a_real_number(laby.get_m());

        secondX = randomX;
        secondY = randomY;

        string where_we_would_do_a_hole = potw[rand() % 4];

        if (where_we_would_do_a_hole == "UP" && randomX != 0) {
            secondX--;
        }
        if (where_we_would_do_a_hole == "DOWN" && randomX != do_a_real_number(laby.get_n())-1){
            secondX++;
        }
        if (where_we_would_do_a_hole == "LEFT" && randomY != 0) {
            secondY--;
        }
        if (where_we_would_do_a_hole == "RIGHT" && randomY != do_a_real_number(laby.get_m())-1) {
            secondY++;
        }
        insert_something_between_two_points_2d(randomX, randomY, secondX, secondY, false);
    }
    
    vector<vector<int>> laby_at_graph() {
        int i, j; int Xsize = laby.get_n(); int Ysize = laby.get_m(); int vert_quantity = do_a_real_number(Xsize)* do_a_real_number(Ysize);
        vector<vector<int>> graph(vert_quantity);

        Xsize = do_a_real_number(Xsize); Ysize = do_a_real_number(Ysize);
        for (i = 0; i < Xsize; i++) {
            for (j = 0; j < Ysize; j++) {
                int ic = do_a_correct_number(i);
                int jc = do_a_correct_number(j);
                bool cur = laby.get_kl(ic, jc);
                if (cur) {
                    if (!(this_is_entry_or_exit(ic - 1, jc)) && laby.get_kl(ic - 1, jc) && laby.get_kl(ic - 2, jc)) {
                            graph[(Ysize * i + j)].push_back((Ysize * (i-1) + j ));
                    }
                    if (!(this_is_entry_or_exit(ic, jc-1)) && laby.get_kl(ic, jc - 1) && laby.get_kl(ic, jc - 2)) {
                            graph[(Ysize * i + j)].push_back((Ysize * i + (j-1) ));
                    }
                    if (!(this_is_entry_or_exit(ic + 1, jc)) && laby.get_kl(ic + 1, jc) && laby.get_kl(ic + 2, jc)) {
                            graph[(Ysize * i + j)].push_back((Ysize * (i+1) + j ));
                    }
                    if (!(this_is_entry_or_exit(ic, jc+1)) && laby.get_kl(ic, jc + 1) && laby.get_kl(ic, jc + 2)) {
                            graph[(Ysize * i + j)].push_back((Ysize * i + (j+1) ));
                    }
                }
            }
        }

        return graph;
    }

    bool this_is_entry_or_exit(int x_coord, int y_coord) {
        return x_coord == start[0] && y_coord == start[1] || x_coord == exit[0] && y_coord == exit[1] ;
    }

    
    void dfs()  {
        int Xsize = laby.get_n(); int Ysize = laby.get_m(); int vert_quantity = do_a_real_number(Ysize) * do_a_real_number(Xsize);
        vector<bool> visited(vert_quantity);
        stack<int> stackan;


        int one_step_start; int one_step_exit;
        vector<int> one_step;

        one_step_start = do_one_step(start[0], start[1]);

        one_step_exit = do_one_step(exit[0], exit[1]);

        int top_of_stackan; vector<int> next;

        stackan.push(one_step_start);

        vector<vector<int>> list_of_edges = laby_at_graph();

        while (!stackan.empty()) {
            top_of_stackan = stackan.top();
            stackan.pop();
            visited[top_of_stackan] = true;
            
            next = list_of_edges[top_of_stackan];
                        
            
            for (int i = 0; i < next.size(); i++) {
                if (!stack_contains_elem(stackan, next[i]) && !visited[next[i]]) stackan.push(next[i]);
            }
        }

        if (visited[one_step_exit]) {
            cout << "Hooray! You find exit"<<endl;
        }
        else {
            cout << "There is no way";
        }
    }
    void bfs() {
        int Xsize = laby.get_n(); int Ysize = laby.get_m(); int vert_quantity = do_a_real_number(Xsize) * do_a_real_number(Ysize);
        vector<bool> visited(vert_quantity);
        deque<int> queue;


        int one_step_start; int one_step_exit;
        vector<int> one_step;

        one_step_start = do_one_step(start[0], start[1]);

        one_step_exit = do_one_step(exit[0], exit[1]);

        int top_of_queue; vector<int> next;


        queue.push_back(one_step_start);

        vector<vector<int>> list_of_edges = laby_at_graph();

        while (!queue.empty()) {
            top_of_queue = queue.front(); 
            queue.pop_front();
            visited[top_of_queue] = true;

            next = list_of_edges[top_of_queue];

            for (int i = 0; i < next.size(); i++) {
                if (!deque_contains_elem(queue, next[i]) && !visited[next[i]]) queue.push_back(next[i]);
            }
        }

        if (visited[one_step_exit]) {
            cout << "Hooray! You find exit"<<endl;
        }
        else {
            cout << "There is no way";
        }
    }
    void output_stack(stack<int> st) {
        while (!st.empty()) {
            st.pop();
        }
    }

    bool deque_contains_elem(deque<int> q, int elem) {
        while (!q.empty()) {
            if (q.front() == elem) return true;
            q.pop_front();
        }
        return false;
    }

    bool stack_contains_elem(stack<int> st, int elem) {
        while (!st.empty()) {
            if (st.top() == elem) return true;
            st.pop();
        }
        return false;
    }

    int do_one_step(int X, int Y) {
        if (X == 0) {
            X++;
        }
        if (X == laby.get_n() - 1) {
            X--;
        }
        if (Y == 0) {
            Y++;
        }
        if (Y == laby.get_m() - 1) {
            Y--;
        }
        return do_a_real_number(laby.get_m()) * do_a_real_number(X) + do_a_real_number(Y);

    }
        

    void deikstra() {
        int Xsize = do_a_real_number(laby.get_n()); int Ysize = do_a_real_number(laby.get_m()); int n = Xsize * Ysize;

        vector<vector<pair<int, int>>> g(n);

        vector<vector<int>> sp_smezh = laby_at_graph();

        for (int i = 0; i < sp_smezh.size(); i++) {
            for (int j = 0; j < sp_smezh[i].size(); j++) {
                g[i].push_back(make_pair(sp_smezh[i][j], 1));
            }
        }

        int s = do_one_step(start[0], start[1]);
        int e = do_one_step(exit[0], exit[1]);


        vector<int> d(n, big_number), p(n);
        d[s] = 0;
        priority_queue<pair<int, int>> q;
        q.push(make_pair(0, s));
        while (!q.empty()) {
            int v = q.top().second, cur_d = -q.top().first;
            q.pop();
            if (cur_d > d[v])  continue;

            for (int j = 0; j < g[v].size(); ++j) {
                int to = g[v][j].first,
                    len = g[v][j].second;
                if (d[v] + len < d[to]) {
                    d[to] = d[v] + len;
                    p[to] = v;
                    q.push(make_pair(-d[to], to));
                }
            }
        }
        if (d[e] != big_number) {
            vector<int> path;
            for (int v = e; v != s; v = p[v])
                path.push_back(v);
            path.push_back(s);

            vector<int> not_sorted_path = path;

            sort(path.begin(), path.end());

            Bool_Matrix with_path = laby;

            for (int i = 0; i < laby.get_n(); i++) {
                for (int j = 0; j < laby.get_m(); j++) {
                    with_path.set_new_value(i, j, 0);
                    int this_coordinate = do_a_real_number(i) * Ysize + do_a_real_number(j);
                    if (i * j % 2 && binary_search(path.begin(),path.end(), this_coordinate) || this_is_entry_or_exit(i,j)) {
                        with_path.set_new_value(i, j, 1);
                    }
                }
            }

            for (int i = 0; i < not_sorted_path.size()-1; i++) {
                int X1coord = not_sorted_path[i] / Ysize;
                int Y1coord = not_sorted_path[i] % Ysize;
                int X2coord = not_sorted_path[i+1] / Ysize;
                int Y2coord = not_sorted_path[i+1] % Ysize;
                with_path.set_new_value((X1coord+X2coord)+1, (Y1coord+Y2coord)+1, 1);
            }

            with_path.output_matrix();

        }
        else {
            cout << "There is no way";
        }
    }

    void wave_alg() {
        int Xsize = laby.get_n(); int Ysize = laby.get_m(); int vert_quantity = do_a_real_number(Xsize) * do_a_real_number(Ysize);
        vector<bool> visited(vert_quantity);
        deque<int> queue;


        int one_step_start; int one_step_exit;
        vector<int> one_step;

        one_step_start = do_one_step(start[0], start[1]);

        one_step_exit = do_one_step(exit[0], exit[1]);

        int top_of_queue; vector<int> next_wave;


        queue.push_back(one_step_start);

        vector<vector<int>> list_of_edges = laby_at_graph();

        while (!queue.empty()) {
            top_of_queue = queue.front(); 
            queue.pop_front();
            visited[top_of_queue] = true;

            next_wave = list_of_edges[top_of_queue];

            for (int i = 0; i < next_wave.size(); i++) {
                if (!deque_contains_elem(queue, next_wave[i]) && !visited[next_wave[i]]) queue.push_back(next_wave[i]);
            }
        }
    }
    

};

int main()
{
    setlocale(LC_ALL, "Russian");

    Labyrinth laby(5,20);
    char c = 0;
    while (c != '1') {
        laby.create_2d_lab();
        laby.create_entry_and_exit();
        laby.output_laby();
        laby.wave_alg();
        laby.dfs();
        laby.bfs(); 
        laby.deikstra();
        cin >> c;
    }
}