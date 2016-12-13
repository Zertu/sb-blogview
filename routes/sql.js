const mysql = require('mysql')
const pool = mysql.createPool({
    host: '104.128.86.50',
    user: 'root',
    password: 'a7864548',
    database: 'sbblog',
    port: 3306,
})
/**
 * 
 * 
 * @param {any} sql字符串
 * @param {any} 回调函数，有err和rows参数
 */
function querysql(sqlstring,fn) {
    pool.getConnection((err,con)=>{
        if(err){
            console.error(err)
        }else{
            con.query(sqlstring,(err,rows)=>{
                console.info(sqlstring)
                if(fn)fn(err,rows)
            })
        }
    })
}
 module.exports={
     querysql:querysql
 }