const express = require('express'),
  router = express.Router(),
  user = require('./user'),
  session = require('express-session')
  /* GET home page. */
let app = express()
app.use(session({
  secret: 'mysbblog', // 建议使用 128 个字符的随机字符串
  cookie: {
    user: {
      user_id: 2,
      userName: 'zertu',
      passWord: 'a7864548',
      Email: '290055513@qq.com',
      lastLoginTime: '2016-12-12T16:00:00.000Z',
      createTime: '1995-06-06T16:00:00.000Z',
      group_Id: 0
    },
    maxAge: 1000 * 60 * 60 
  }
}))
router.get('/', function (req, res, next) {
  res.render('index', {
    title: 'sb-blog'
  });
});
router.get('/article', (req, res) => {
  res.render('article', {
    title: 'sb-blog',
    session: req.session
  })
})


router.post('/login', (req, res) => {
  let row = user.checklogin(req, res, function (rows) {
    res.render('article', {
      title: 'sb-blog',
      session: rows
    })
  })

})

module.exports = router;